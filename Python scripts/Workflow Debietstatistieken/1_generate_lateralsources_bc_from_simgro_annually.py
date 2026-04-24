"""
Generate DFM_lateral_sources.bc from SIMGRO output (sw_dtsw or sw_dtgw)
Base script: "1_1_generate_lateralsources_bc_from_simgro_output.py"
Location: reken05 - "D:\Bovenregionale_stresstest_ARKNZK\0_scripts"

Modification: annual files

Author:         Petra Hulsman (HDSR)
Last update:    13/04/2026
Running environment: modflow01 - env_lats_Petra
"""


import pandas as pd
import numpy as np
import os
from datetime import datetime as dt
from datetime import timedelta as td
from tqdm.auto import tqdm

class paths():
    root        = r'D:/Petra/UGM2lat'
    simgro_path = r'D:\hydromedah_P23006_KALIB_V4B_LONGPERIOD_20241024_20110101_20240101_250m\output\metaswap'
    path_orig   = root + r'/Data_input/DFM_lateral_sources_orig.bc'
    dir_interm  = root + r'/Data_intermediate'
    
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
        val = float(np.fromfile(f,dtype=np.float32,count=1))
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


def read_bda(simgro_path, filename, bytesize=8, year_filter=None, locind=None):     
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

def replace_in_orig_bc_file(path_new,y,season,tref):
    ''' take the original bc file, replace the newly generated data from this script, and keep data for remaining locations'''
    # Define additional paths
    path_repl = paths.root + r'/Data_output/DFM_lateral_sources_' + season + str(y) + '.bc'
    path_ext  = paths.root + r'/Data_intermediate/DFM_lateral_sources_extra_' + season + str(y) + '.bc'
    
    # Load data
    lats_orig = open(paths.path_orig, encoding='utf-8', errors='ignore')
    lats_new  = open(path_new, encoding='utf-8', errors='ignore')
    lines_orig= lats_orig.readlines()
    lines_new = lats_new.readlines()
    
    # create new output file
    if os.path.exists(path_repl): os.remove(path_repl)
    if os.path.exists(path_ext): os.remove(path_ext)
    df_out = open(path_repl, 'w')
    df_ext = open(path_ext, 'w')
    
    # check for groupname types
    check1 = np.unique(np.array([l for l in lines_orig if l.startswith('[')]))
    check2 = np.unique(np.array([l for l in lines_new if l.startswith('[')]))
    
    # start index blocks of info
    groupname       = '[forcing]\n'
    idx_blocks_new  = np.array([i for i in tqdm(range(0,len(lines_new))) if lines_new[i]==groupname])
    names_new       = np.array([lines_new[i+1].split('=')[1].replace('\n','').strip() for i in idx_blocks_new])
    idx_blocks_orig = np.array([i for i in tqdm(range(0,len(lines_orig))) if lines_orig[i]==groupname])
    names_orig      = np.array([lines_orig[i+1].split('=')[1].replace('\n','').strip() for i in idx_blocks_orig])
    
    # keep first block
    for ii in range(0,idx_blocks_orig[0]):
        df_out.write(lines_orig[ii])
    
    with tqdm(total=len(idx_blocks_orig)-1) as pbar:
        for i in range(0,len(idx_blocks_orig)):
        # for i in range(1286,len(idx_blocks_orig)):
        # for i in range(0,10):
            if names_orig[i] not in names_new:    
                # use orig data
                row_start = int(idx_blocks_orig[i])
                if row_start != idx_blocks_orig[-1]:
                    row_end   = idx_blocks_orig[i + 1 ]
                else:
                    row_end = len(lines_orig)                              
                
                for ii in range(row_start,row_end):
                    # change dates
                    lines_orig[row_start:row_end]
                    if 'since' in lines_orig[ii]:
                        lines_orig[ii] = lines_orig[ii].replace('2022-01-01 00:00:00',tref)
                        lines_orig[ii+4] = lines_orig[ii+4].replace('525540','31536000')                        
                    df_ext.write(lines_orig[ii])
                    df_out.write(lines_orig[ii])
            else:
                # use new data
                ii = int(np.where(np.array(names_new)==names_orig[i])[0][0])
                row_start = idx_blocks_new[ii]
                if row_start != idx_blocks_new[-1]:
                    row_end   = idx_blocks_new[ii + 1 ] 
                else:
                    row_end = len(lines_new)                              
                
                for iii in range(row_start,row_end): 
                    sstring = 'seconds since'
                    
                    if sstring in lines_new[iii]:
                        ref_year = lines_new[iii].split('=')[1].replace('seconds since ','')[0:5]
                        ref_year = ref_year.strip()
                        if ref_year!=tref[0:4]:
                            print('Error, i=', i, 'year = ', ref_year)
                            # error
                        lines_new[iii] = lines_new[iii]
                        
                    df_out.write(lines_new[iii])
            pbar.update(1)
    
    # close files
    lats_orig.close()
    lats_new.close()
    df_out.close()
    df_ext.close()
    
    
if __name__ == "__main__":
    print('Start')
    
    # load laterals from SIMGRO output
    msw_file = 'sw_dtsw'
    df = read_bda(
        simgro_path=paths.simgro_path, 
        filename=msw_file
    )
    
    # do some preprocessing and save to csv file    
    lateral_format_string = 'AW{lateral:04.0f}'
    df.columns = [lateral_format_string.format(lateral=int(c)) for c in df.columns]
    df.dropna(how='all', axis=1, inplace=True)
    df.to_csv(os.path.join(paths.dir_interm, 'DFM_lateral_sources.csv'), sep=';', index_label='datetime')
    
    
    # df = pd.read_csv(os.path.join(output_path, 'DFM_lateral_sources.csv'), delimiter=';').set_index('datetime')
    # df.index = pd.to_datetime(df.index)
    
    # split periods
    for y in range(2014,2024):
        print('Splitting by season', y)
        
        # split by season
        # winter: 15 okt t/m 14 apr
        tw_1 = str(y) + '-10-15 00:00:00'
        tw_2 = str(y+1) + '-04-15 00:00:00'
        # zomer: 15 apr t/m 14 okt
        tz_1 = str(y) + '-04-15 00:00:00'
        tz_2 = str(y) + '-10-15 00:00:00'
        
        # select dataframe
        df_z      = df[(df.index >= tz_1) & (df.index <= tz_2)]
        df_w      = df[(df.index >= tw_1) & (df.index <= tw_2)]
        
        # define output paths        
        path_new_z = os.path.join(paths.dir_interm, 'DFM_lateral_sources_UGM_zomer'+str(y)+'.bc')
        path_new_w = os.path.join(paths.dir_interm, 'DFM_lateral_sources_UGM_winter'+str(y)+'.bc')
        
        # write lateral sources BC file
        write_lateral_sources_bc_file(laterals_df=df_z,bc_file=path_new_z)
        write_lateral_sources_bc_file(laterals_df=df_w,bc_file=path_new_w)
        
        # replace in orig bc file
        replace_in_orig_bc_file(path_new_z,y,'zomer',tz_1)
        replace_in_orig_bc_file(path_new_w,y,'winter',tw_1)
        
print('Done')
