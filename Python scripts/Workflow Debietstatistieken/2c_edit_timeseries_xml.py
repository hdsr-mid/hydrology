# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 26/02/2026
virtual environment used: geo-env

"""

import numpy as np
import os
from tqdm import tqdm
import pandas as pd

class paths():
    path_in = r'D:/workingdir/1_InputData/ModelData/timeseries_import.txt'
    path_out= r'D:/workingdir/3_Output/ModelData/timeseries_import'
    

if __name__ == "__main__":
    
    # Load data
    data_in = open(paths.path_in, encoding='utf-8', errors='ignore')
    lines   = data_in.readlines()
    
    # create new output file
    path_out_w = paths.path_out + '_winter.xml'
    path_out_z = paths.path_out + '_zomer.xml'
    if os.path.exists(path_out_w): os.remove(path_out_w)
    if os.path.exists(path_out_z): os.remove(path_out_z)
    data_out_w = open(path_out_w, 'w')
    data_out_z = open(path_out_z, 'w')
    
    # start index blocks of info
    groupname   = '  <series>\n'
    idx_blocks  = np.array([i for i in tqdm(range(0,len(lines))) if lines[i]==groupname])
    
    # keep first block
    for ii in range(0,idx_blocks[0]):
        data_out_w.write(lines[ii])    
        data_out_z.write(lines[ii])    
    
    with tqdm(total=len(idx_blocks)-1) as pbar:
        for ii in range(0,len(idx_blocks)):
        # for ii in range(0,5):
            # start/end row index
            row_start = idx_blocks[ii]
            if row_start != idx_blocks[-1]:
                row_end   = idx_blocks[ii + 1 ] 
            else:
                row_end = len(lines)                              
            
            # get winter/summer data
            groupname   = '    <event '
            idx_data    = np.array([i for i in range(row_start,row_end) if lines[i].startswith(groupname)])
            dates       = [lines[i].split('"')[1] for i in idx_data]
            times       = [lines[i].split('"')[3] for i in idx_data]
            values      = [lines[i].split('"')[5] for i in idx_data]
            dates       = [pd.to_datetime(dates[t]+' '+times[t]) for t in range(0,len(dates))]
            df_i        = pd.DataFrame({'Dates':dates,'value':values})
            if len(np.unique(np.array(values)))==1:
                value_w = 'value="' + str(float(np.unique(np.array(values)))) + '"'
                value_z = 'value="' + str(float(np.unique(np.array(values)))) + '"'
            else:
                cond = np.array([pd.to_datetime(d).year for d in df_i.Dates.values])==2022
                df_i = df_i[cond]
                mask_z      = (df_i['Dates'] >= '2022-04-15 00:00:00') & (df_i['Dates'] <= '2022-10-14 23:00:00')
                mask_w      = (df_i['Dates'] <= '2022-04-14 23:00:00') | (df_i['Dates'] >= '2022-10-15 00:00:00')
                value_w     = 'value="' + str(float(np.unique(df_i.loc[mask_w].value))) + '"'
                value_z     = 'value="' + str(float(np.unique(df_i.loc[mask_z].value))) + '"'
                
            # save lines
            basis_line    = lines[idx_data[0]]
            value2replace = 'value="'+basis_line.split('"')[5]+'"'
            for jj in range(row_start,row_end):
                if jj == idx_data[0]:
                    data_out_w.write(basis_line.replace('2022','2011').replace(value2replace,value_w))
                    data_out_z.write(basis_line.replace('2022','2011').replace(value2replace,value_z))
                if jj == idx_data[-1]:
                    data_out_w.write(basis_line.replace('2022','2025').replace(value2replace,value_w))
                    data_out_z.write(basis_line.replace('2022','2025').replace(value2replace,value_z))
                if 'startDate' in lines[jj]:
                    lines[jj] = lines[jj].replace('2022-01-01','2011-01-01')
                    # some manual exceptions:
                    lines[jj] = lines[jj].replace('2022-01-26','2011-01-01')
                        
                if 'endDate' in lines[jj]:
                    lines[jj] = lines[jj].replace('2022-12-31','2025-01-01')
                    # some manual exceptions:
                    lines[jj] = lines[jj].replace('2022-03-01','2025-01-01')
                if jj not in idx_data:
                    data_out_w.write(lines[jj])
                    data_out_z.write(lines[jj])
                
            pbar.update(1)
    
    # close files
    data_in.close()
    data_out_w.close()
    data_out_z.close()
    
    
    