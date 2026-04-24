# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 19/02/2026
virtual environment used: geo-env

"""

import numpy as np
import os
from tqdm import tqdm
import pandas as pd

class paths():
    path_in = r'D:/workingdir/1_InputData/ModelData//FlowFM_structures.bc'
    path_out= r'D:/workingdir/3_Output/ModelData/FlowFM_structures'
    

if __name__ == "__main__":
    
    # Load data
    data_in = open(paths.path_in, encoding='utf-8', errors='ignore')
    lines   = data_in.readlines()
    
    # create new output file
    path_out_w = paths.path_out + '_winter.bc'
    path_out_z = paths.path_out + '_zomer.bc'
    if os.path.exists(path_out_w): os.remove(path_out_w)
    if os.path.exists(path_out_z): os.remove(path_out_z)
    data_out_w = open(path_out_w, 'w')
    data_out_z = open(path_out_z, 'w')
    
    # start index blocks of info
    groupname   = '[forcing]\n'
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
            idx_data    = np.array([i for i in range(row_start,row_end) if '=' not in lines[i] and groupname not in lines[i] and lines[i]!='\n'])
            idx_tref    = int(np.array([i for i in range(row_start,row_end) if 'since' in lines[i]]))
            dates       = [int(lines[i].strip().split(' ')[0]) for i in idx_data]
            values      = [float(lines[i].strip().split(' ')[1]) for i in idx_data]
            tref        = lines[idx_tref].strip().split('=')[1].replace(' minutes since ','').replace('\n','')
            tref        = pd.to_datetime(tref)
            if 'minutes' in lines[idx_tref]:
                dates       = [tref+pd.Timedelta(minutes=t) for t in dates]
            else:
                dates = ''
                print('Error time unit != minutes')
                break
                
            if (tref.year!=2022) & (tref.year!=2011):
                print('Error ref year: ', tref.year, ' =! 2022/2011')
                break
                
            df_i        = pd.DataFrame({'Dates':dates,'value':values})
            if len(np.unique(np.array(values)))==1:
                value_w = str(float(np.unique(np.array(values))))
                value_z = str(float(np.unique(np.array(values))))
            else:
                cond = np.array([pd.to_datetime(d).year for d in df_i.Dates.values])==tref.year
                df_i = df_i[cond]
                mask_z      = (df_i['Dates'] >= str(tref.year)+'-04-15 00:00:00') & (df_i['Dates'] <= str(tref.year)+'-10-14 23:00:00')
                mask_w      = (df_i['Dates'] <= str(tref.year)+'-04-14 23:00:00') | (df_i['Dates'] >= str(tref.year)+'-10-15 00:00:00')
                value_w     = str(float(np.unique(df_i.loc[mask_w].value)))
                value_z     = str(float(np.unique(df_i.loc[mask_z].value)))
                
            # new lines
            lines[idx_tref]= lines[idx_tref].replace('2022','2011')
            new_lines_w1   = '    0 ' + str(value_w) + ' \n'
            new_lines_w2   = '    410227200 ' + str(value_w) + ' \n'
            new_lines_z1   = '    0 ' + str(value_z) + ' \n'
            new_lines_z2   = '    410227200 ' + str(value_z) + ' \n'
            
            # save lines
            for jj in range(row_start,row_end):
                if jj == idx_data[0]:
                    data_out_w.write(new_lines_w1)
                    data_out_z.write(new_lines_z1)
                if jj == idx_data[-1]:
                    data_out_w.write(new_lines_w2)
                    data_out_z.write(new_lines_z2)                    
                if jj not in idx_data:
                    data_out_w.write(lines[jj])
                    data_out_z.write(lines[jj])
                
            pbar.update(1)
    
    # close files
    data_in.close()
    data_out_w.close()
    data_out_z.close()
    
    
    