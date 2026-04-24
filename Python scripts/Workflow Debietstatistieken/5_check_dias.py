# -*- coding: utf-8 -*-
"""
@author     : MaaS-user, Petra Hulsman
Last update : 13/04/2026
virtual environment used: geo-env

"""

import numpy as np
import os
from tqdm import tqdm
import pandas as pd

class paths():
    root    = r'D:/workingdir/Model_runs/dia'
    path_out1  = r'D:/workingdir/3_Output/txt/DFM_knijpen_temp.txt'
    path_out2  = r'D:/workingdir/3_Output/txt/DFM_changes_temp.txt'
    path_out_dt = r'D:/workingdir/3_Output/txt/output_tijdstap.txt'
    path_out_kn = r'D:/workingdir/3_Output/txt/output_knijpen.txt'
    path_out_ch = r'D:/workingdir/3_Output/txt/output_changes.txt'
    path_out_ch2= r'D:/workingdir/3_Output/txt/output_changes_samenvatting.txt'
    
    

if __name__ == "__main__":
    
    files = [f for f in os.listdir(paths.root) if '.dia' in f]
    
    # create new output file
    df_dt = pd.DataFrame()
    df_kn = pd.DataFrame()
    df_ch = pd.DataFrame()
    
    data_out_ch_temp  = open(paths.path_out2, 'w')
    
    for file in files:
        print(file)
        
        # Load data
        data_in    = open(os.path.join(paths.root,file), encoding='utf-8', errors='ignore')
        lines      = data_in.readlines()
        data_out1  = open(paths.path_out1, 'w')
        
        # settings
        n_knijpen = 0
        n_changes = 0
        dt_list   = []
        
        for ii in tqdm(range(0,len(lines))):
        # for ii in range(0,5):
            if ('WARNING' in lines[ii]) & ('Discharge' in lines[ii]) & ('limited' in lines[ii]):
                data_out1.write(lines[ii])
                n_knijpen = n_knijpen + 1
            if ('WARNING' in lines[ii]) & ('changed' in lines[ii]):
                data_out_ch_temp.write(lines[ii])                
            if ('INFO' in lines[ii]) & ('%' in lines[ii]) & ('Sim. time' not in lines[ii]) & ('?' not in lines[ii]):
                line_parts = [l.strip() for l in lines[ii].split('   ')]
                percentage = float(line_parts[-2].replace('%',''))
                if percentage>0:
                    dt_list.append(float(line_parts[-1]))
        
        # close files
        data_in.close()
        data_out1.close()
        del lines, data_in, data_out1
        
        # file details
        Year      = int(file.replace('DFM_','').replace('winter','').replace('zomer','').replace('.dia','').replace('stopped',''))
        Season    = file.replace('DFM_','').replace('.dia','').replace('stopped','').replace(str(Year),'')
        
        
        # summary: knijpen
        remove = ['** WARNING: Discharge through pump ',
                  ' is limited below capacity by water volume on suction side.\n']
        data_in   = open(paths.path_out1, encoding='utf-8', errors='ignore')
        lines     = data_in.readlines()
        lines_u   = np.unique(lines)
        l_count   = [len(np.where(np.array(lines)==l)[0]) for l in tqdm(lines_u)]
        
        for ii in tqdm(range(0,len(lines_u))): 
            line_ii = lines_u[ii]
            for r in remove: line_ii = line_ii.replace(r,'')
            perc = l_count[ii]/n_knijpen    
            df_kn = pd.concat([df_kn,pd.DataFrame({'Year':[Year],'Season': [Season],'Location':[line_ii],'Occurences':[perc]})],ignore_index=True)            
        data_in.close()
        os.remove(paths.path_out1)
                
        # summary: timestep
        trange = np.arange(0,70,5)
        for idx in range(1,len(trange)):
            perc = len(np.where((np.array(dt_list) >= trange[idx-1]) & (np.array(dt_list) <= trange[idx]))[0])/len(dt_list)
            df_dt = pd.concat([df_dt,pd.DataFrame({'Year':[Year],'Season': [Season],'Timestep':[' ' + str(trange[idx-1]) + ' - ' + str(trange[idx])],'Occurences':[perc]})],ignore_index=True)            
                
        
    # save dataframe
    df_dt.to_csv(paths.path_out_dt)
    df_kn.to_csv(paths.path_out_kn)
    
    
    # summary: changes
    data_out_ch_temp.close()
    data_in   = open(paths.path_out2, encoding='utf-8', errors='ignore')
    lines     = data_in.readlines()
    lines_u   = np.unique(lines)
    data_out  = open(paths.path_out_ch, 'w')
    for ii in tqdm(range(0,len(lines_u))): data_out.write(lines_u[ii])
    data_in.close()
    os.remove(paths.path_out2)
    