# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 10:09:53 2025

@author: MaaS-user

Stappen:
    1. identificeer blokken in ini bestand
    2. vind id_name
    3. behoud blok als id_name terug te vinden is in het netcdf bestand (network_branch_id)
    uitzondering voor crsdef:
    3. behoud blok als id_name terug te vinden is in het crsloc.ini bestand (Id)
    
"""

import numpy as np
import pandas as pd
import os
from tqdm import tqdm
import sys
import xarray as xr

    
def clip_ini(path,fn_in, groupname):
    print('clip ' + fn_in)
    # load data
    df_in  = open(os.path.join(path,fn_in))
    lines  = df_in.readlines()
    
    # create new output file
    fn_out   = fn_in.replace('.ini','_sel.ini')
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    # selection based on netcdf
    xds       = xr.open_dataset(os.path.join(path,'network.nc'))
    branch_id = xds.network_branch_id.values.astype(str)
    branch_id = np.array([x.replace(' ','') for x in branch_id ])
    # print(branch_id)
    
    # uitzondering: selecteer o.b.v. crsloc
    if os.path.basename(fn_in) == 'crsdef.ini':
        # based on crsloc
        df_loc = open(os.path.join(path,'crsloc.ini'))
        l_loc  = df_loc.readlines()
        IDs1   = [l_loc[x].split('=')[1] for x in np.where(np.array(l_loc)=='[CrossSection]\n')[0] + 5]
        IDs1   = np.array([x.replace(' ','')[:-1] for x in IDs1 ])
        df_loc.close()
        
        # based on structures
        df_loc = open(os.path.join(path,'structures.ini'))
        l_loc  = df_loc.readlines()
        IDs2   = [l_loc[x].split('=')[1] for x in np.where(np.array(l_loc)=='[Structure]\n')[0] + 1]
        IDs2   = np.array([x.replace(' ','')[:-1] for x in IDs2 ])
        df_loc.close()
        
        # combine
        IDs = np.append(IDs1,IDs2)
        
    # start index blocks of info
    idx_blocks = np.where(np.array(lines)==groupname)[0] 
    
    # check for groupname types
    check = np.unique(np.array([l for l in lines if '=' not in l]))
    
    # keep first block
    for ii in range(0,idx_blocks[0]):
        df_out.write(lines[ii])
    
    with tqdm(total=len(idx_blocks)-1,desc=fn_in) as pbar:
        for i in range(0,len(idx_blocks)):
            row_start = idx_blocks[i]
            if row_start != idx_blocks[-1]:
                row_end   = idx_blocks[i + 1 ] -1 
            else:
                row_end = len(lines)
            
            # keep next block under specific conditions
            categories = np.array([l.split('=')[0].replace(' ','') for l in lines[row_start:row_end]])
            if os.path.basename(fn_in) == 'crsdef.ini':
                idx_crit   = row_start + np.where(categories=='Id')[0][0] # row in block starting with "Id"
                id_name    = lines[idx_crit].split('=')[1].replace(' ','')[:-1] # corresponding string
                cond       = id_name in IDs                
            else:
                idx_crit   = row_start + np.where(categories=='branchId')[0][0] # row in block starting with "branchID"
                id_name    = lines[idx_crit].split('=')[1].replace(' ','')[:-1] # corresponding string
                cond       = id_name in branch_id

            if cond:                  
                for ii in range(row_start,row_end):
                    df_out.write(lines[ii])
                
            pbar.update(1)
    
    df_in.close()
    df_out.close()

    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
    

    
