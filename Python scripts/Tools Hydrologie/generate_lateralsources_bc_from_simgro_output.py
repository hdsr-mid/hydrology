"""
Generate DFM_lateral_sources.bc from SIMGRO output (sw_dtsw or sw_dtgw)
Copied functions read_key (ReadKey), get_timeseries (GetTimeSeries_V2), read_bda (Read_BDA_V2) 
from "Hydrologische Informatieproducten" toolbox

Author:         Boyan Domhof (HDSR)
Last update:    29/11/2024
"""


import pandas as pd
import numpy as np
import os
from datetime import datetime as dt
from datetime import timedelta as td
from tqdm.auto import tqdm


def read_key(filename):
        """ Functie om een SIMGRO key file in te lezen"""
        kf = open(filename+'.key','r')
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


def get_timeseries(f, numtim, numlocs, numvar, loc, var, bytesize, timeindexfilter=None):
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
        val = float(np.fromfile(f,dtype=np.float32,count=1))
        timeserie[tim] = val
    return timeserie  


def read_tim(filename):
    """ Functie om een SIMGRO tim file in te lezen"""
    tf = open(filename+'.tim','r')
    lines = tf.readlines()
    tf.close()
    times = [(float(line.lstrip().split(' ')[0]), int(line.lstrip().split(' ')[2]))  for line in lines]
    dates = [dt.strptime(str(time[1])+'/01/01 00:00:00','%Y/%m/%d %H:%M:%S') + td(days=time[0]) for time in times]
    dates = [date+td(0,1) if date.minute==59 else date for date in dates]
    return dates


def read_bda(simgro_path, filename, bytesize=None, year_filter=None, locind=None):     
    """ Functie om een binair SIMGRO bestand uit te lezen"""
    # To calculate lateral: -[Vdsreg] - [Vsues]
    variables = ('Vdsreg','Vsues')

    (locations, variabelen, period, key_bytesize) =  read_key(os.path.join(simgro_path, filename))        
    (timesteps) = read_tim(os.path.join(simgro_path, filename))
    if bytesize is None:
        bytesize = key_bytesize
    tijdstap = np.floor((timesteps[1]-timesteps[0]).total_seconds())
    if tijdstap == 3600.:
        frequency = 'H'
    if tijdstap == 86400.:
        frequency = 'D'
    varind = [variabelen.index(var) for var in variables]
    
    # create neat time series from, with correct interval and proper length
    # takes into accoutn period type (as specified in the key file)
    if period == 1:
        timesteps = [timesteps[0] + td(seconds=tijdstap*x) for x in range(len(timesteps)-1)]
    else:
        timesteps = [timesteps[0] + td(seconds=tijdstap*x) for x in range(len(timesteps))]
    timesteps = list(timesteps)
    
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
            f = open(os.path.join(simgro_path,filename+'.bda'), "rb")
            ts.append(get_timeseries(f, numtim, numpts, numvar, locind, vi, bytesize=bytesize, timeindexfilter=timefilter_index))
            f.close()
        # "(-ts[0][i] - ts[1][i])" equals "-[Vdsreg] - [Vsues]"
        df.loc[:,location] = [(-ts[0][i] - ts[1][i])/tijdstap for i in range(len(ts[0]))]
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
            f.write(header.replace('$LATERAL$', c).replace('$TSTART$', f'{laterals_df.index[0]:%Y-%m-%d %H:%M:%S}'))
            for t, v in zip(laterals_df.index, laterals_df[c].values):
                f.write(f'    {int((t-laterals_df.index[0]).total_seconds())}    {v:.6f}\n')


print('Start')

# load laterals from SIMGRO output
simgro_path = r'd:\Bovenregionale_stresstest_ARKNZK\1_basis_scenario\hydromedah_output'
msw_file = 'sw_dtsw'
df = read_bda(
    simgro_path=simgro_path, 
    filename=msw_file, 
    bytesize=8
)

# do some preprocessing and save to csv file
lateral_format_string = 'AW{lateral:04.0f}'
df.columns = [lateral_format_string.format(lateral=int(c)) for c in df.columns]
df.dropna(how='all', axis=1, inplace=True)
df.to_csv(os.path.join(simgro_path, 'DFM_lateral_sources.csv'), sep=';', index_label='datetime')

# write lateral sources BC file
write_lateral_sources_bc_file(
    laterals_df=df,
    bc_file=os.path.join(simgro_path, 'DFM_lateral_sources.bc')
)

print('Done')
