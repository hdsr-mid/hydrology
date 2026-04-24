"""
Generate DFM_lateral_sources.bc from SIMGRO output (sw_dtsw or sw_dtgw)
Base script: "1_1_generate_lateralsources_bc_from_simgro_output.py"
Location: reken05 - "D:\Bovenregionale_stresstest_ARKNZK\0_scripts"

Author:         Petra Hulsman (HDSR)
Last update:    09/04/2026
"""


import pandas as pd
import numpy as np
import os
from datetime import datetime as dt
from datetime import timedelta as td
from tqdm.auto import tqdm


def read_key(filename):
        """ Functie om een SIMGRO key file in te lezen"""
        kf = open(filename + '.key', 'r')
        lines = kf.readlines()
        kf.close()
        for iline,line in enumerate(lines):
            if line.startswith('*'): continue
            if line.startswith('FORMAT'):
                bytesize = int(line.split('=')[1].strip()[1])
            if line.startswith('PERIOD'):
                period = int(line.split('=')[1].strip())
            if line.startswith('NUMVAR'):
                numvar = int(line.split('=')[1].strip())        
                varlines = lines[iline+3:iline+3+numvar]
                variables = [var.split(' ')[0].strip() for var in varlines]            
            if line.startswith('NUMPTS'):
                numlocs = int(line.split('=')[1].strip())
                loclines = lines[iline+3:iline+3+numlocs]
                locations = [loc.lstrip().split(' ')[0] for loc in loclines]
        return locations, variables, period, bytesize


def get_timeseries(f, numtim, numlocs, numvar, loc, var, bytesize=8, timeindexfilter=None):
    """ Functie om een tijdserie op te halen uit een binair SIMGRO bestand"""
    byte_index = range(var * bytesize, numtim * numlocs * numvar * bytesize, numlocs * numvar * bytesize)
    byte_index = [x + (loc * bytesize * numvar) for x in byte_index]
    # apply timeseries filter
    if timeindexfilter is not None:
        byte_index_filtered = [byte_index[x] for x in timeindexfilter]
        byte_index = list(byte_index_filtered)
    # get the values at the byte indices
    timeserie = [None] * len(byte_index)
    for tim in range(len(byte_index)):
        f.seek(byte_index[tim],0)
        val = float(np.fromfile(f,dtype=np.float32,count=1)[0])
        timeserie[tim] = val
    return timeserie


def read_tim(filename):
    """ Functie om een SIMGRO tim file in te lezen"""
    tf = open(filename + '.tim', 'r')
    lines = tf.readlines()
    tf.close()
    times = [(float(line.lstrip().split(' ')[0]), int(line.lstrip().split(' ')[2]))  for line in lines]
    dates = [dt.strptime(str(time[1]) + '/01/01 00:00:00', '%Y/%m/%d %H:%M:%S') + td(days=time[0]) for time in times]
    dates = [date + td(seconds=1) if ((date.minute == 59) and (date.second == 59)) else date for date in dates]
    # automatically correct inconsistencies in dates
    dates = [dt(year=d.year, month=d.month, day=d.day, hour=d.hour) for d in dates]  # round dates to hour
    return dates


def read_bda(simgro_path, filename, sel_loc, bytesize=8, year_filter=None, locind=None):     
    """ Functie om een binair SIMGRO bestand uit te lezen"""
    # To calculate lateral: -[Vdsreg] - [Vsues]
    variables = ('Vdsreg','Vsues')

    (locations, variabelen, period, key_bytesize) = read_key(os.path.join(simgro_path, filename))        
    (timesteps) = read_tim(os.path.join(simgro_path, filename))
    #if bytesize is None:
    #    bytesize = key_bytesize
    tijdstap = np.floor((timesteps[1]-timesteps[0]).total_seconds())
    if tijdstap == 3600.:
        frequency = 'h'
    if tijdstap == 86400.:
        frequency = 'd'
    varind = [variabelen.index(var) for var in variables]
    
    # create neat time series from, with correct interval and proper length
    # takes into accoutn dperiod type (as specified in the key file)
    if period == 1:
        timesteps = [timesteps[0] + td(seconds=tijdstap*x) for x in range(len(timesteps)-1)]
    else:
        timesteps = [timesteps[0] + td(seconds=tijdstap*x) for x in range(len(timesteps))]
    timesteps = list(timesteps)
    
    # select location
    locations = [l for l in locations if l in sel_loc]
    numpts = len(locations)
    
    # get number of locations, variables and timesteps
    numpts = len(locations)
    numvar = len(variabelen)
    numtim = len(timesteps)

    # apply year filter
    if year_filter is not None:
        timesteps_years = [dt.strftime(tim, '%Y') for tim in timesteps]
        timefilter_index = [i for i, x in enumerate(timesteps_years) if x == year_filter]
        timesteps = list(map(timesteps.__getitem__, timefilter_index))
    else:
        timefilter_index = None
        
    # loop over all locations
    df = pd.DataFrame(np.zeros((numtim, numpts))*np.nan)
    df.columns = locations
    df.index = pd.date_range(timesteps[0], timesteps[-1], freq=frequency)
    
    # get time series for a specific location (index)
    if locind is not None:
        locations = [locind]
    
    
    for location in tqdm(locations, total=numpts):
        if location == '0':
            # skip swnr 0
            continue
        locind = locations.index(location)
        ts = []
        for vi in varind:
            f = open(os.path.join(simgro_path,filename + '.bda'), "rb")
            ts.append(get_timeseries(f, numtim, numpts, numvar, locind, vi, bytesize=bytesize, timeindexfilter=timefilter_index))
            f.close()
        # "(-ts[0][i] - ts[1][i])" equals "-[Vdsreg] - [Vsues]"
        df.loc[:,location] = [(-ts[0][i] - ts[1][i])/tijdstap for i in range(numtim)]
        
    return(df)


def write_lateral_sources_bc_file(
        laterals_df: pd.DataFrame, 
        bc_file: str = 'DFM_lateral_sources.bc',
):
    with open(bc_file, 'w') as f:
        # write general block
        f.write("""[General]
    fileVersion           = 1.01                
    fileType              = boundConds
""")
        # write timeseries blocks
        header = """
[forcing]
    name                  = $LATERAL$              
    function              = timeseries          
    timeInterpolation     = linear              
    quantity              = time                
    unit                  = seconds since $TSTART$
    quantity              = lateral_discharge   
    unit                  = m³/s
""" 
        for c in laterals_df.columns:
            f.write(header.replace('$LATERAL$', c).replace('$TSTART$', f'{laterals_df.index[0]:%Y-%m-%d %H:%M:%S}').replace('2011-', '2021-'))
            for t, v in zip(laterals_df.index, laterals_df[c].values):
                f.write(f'    {int((t-laterals_df.index[0]).total_seconds())}    {v:.6f}\n')



