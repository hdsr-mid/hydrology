# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 11:43:30 2025

@author: MaaS-user
"""


import os
import numpy as np
import xarray as xr
from tqdm import tqdm

def round_netcdf(path):
    print('round netcdf')
    
    # load data
    var_dat = 'Mesh2d_node_z'
    xds     = xr.open_dataset(os.path.join(path,'network.nc'))
    
    # round
    xds[var_dat] = np.round(xds[var_dat],2)
    
    # save file
    xds.to_netcdf(os.path.join(path,'network_2.nc'))
    
    # remove old file, rename new
    xds.close()    
    os.remove(os.path.join(path,'network.nc'))
    os.rename(os.path.join(path,'network_2.nc'),os.path.join(path,'network.nc'))

def round_pliz(path):
    print('round pliz')
    
    # load data
    fn_in   = 'objects.pli_fxw.pliz'
    path_in = os.path.join(path,fn_in)
    df_in   = open(path_in)
    lines   = df_in.readlines()
    
    # create new output file
    fn_out   = fn_in.replace('.pliz','_2.pliz')
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    ll = np.unique([len(line.split(' ')) for line in lines])
    
    with tqdm(total=len(lines)-1,desc=fn_in) as pbar:
        for line in lines:
            if len(line.split(' '))==17:
                line_split = line.split(' ')
                for i in np.arange(4,len(line_split),2):
                    line_split[i] = str(np.round(float(line_split[i]),2))
                line_new   = ' '.join(line_split) + '\n'
                
                df_out.write(line_new)
            else:
                df_out.write(line)
            pbar.update(1)
            
    
    df_in.close()
    df_out.close()

    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))

        
        