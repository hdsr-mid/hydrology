# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 05/05/2026
virtual environment used: geo-env
"""

import os
import xarray as xr
import warnings
warnings.filterwarnings("ignore")


class general():
    model_dir           = r'D:/workingdir/2_Model_runs/netcdf'
    outdir              = r'D:/workingdir/3_Output'
    
    path_base = os.path.join(outdir,'netcdf','base.nc')
    path_WL   = os.path.join(outdir,'netcdf','WL.nc')
    path_Q    = os.path.join(outdir,'netcdf','Q.nc')
    path_all_WL  = os.path.join(outdir,'netcdf','DFM_all_Waterstand.nc')
    path_all_Q   = os.path.join(outdir,'netcdf','DFM_all_Debiet.nc')

def get_base_data():
    path      = os.path.join(general.model_dir,'DFM_map_zomer2014.nc')
    xds       = xr.open_dataset(path)
    var_sel   = [v for v in list(xds.head()) if 'time' not in xds[v].dims and 'projected' not in v and 'mesh1d_edge_branch' not in v]
    xds[var_sel].to_netcdf(general.path_base)
    xds.close()

def get_WL_data():
    seasons = ['zomer','winter']
    for y in range(2014,2024):
        for season in seasons:
            model     = season + str(y)
            if season == 'zomer':
                period  = slice(str(y) + '-04-15',str(y) + '-10-14')
            else:
                period  = slice(str(y) + '-10-15',str(y+1) + '-04-14')
            print(model)
            path      = os.path.join(general.model_dir,'DFM_map_' + season + str(y) + '.nc')
            xds_y     = xr.open_dataset(path).sel(time=period)[['mesh1d_s1']]
            if season + str(y) != 'winter2023': xds_y     = xds_y.sel(time=period)
            if model == 'zomer2014':
                xds_var   = xds_y
            else:
                xds_var   = xr.concat([xds_var,xds_y],dim='time')
            xds_y.close()    
            del xds_y
    xds_var.to_netcdf(general.path_WL)
    xds_var.close()
    

def get_Q_data():
    seasons = ['zomer','winter']
    for y in range(2014,2024):
        for season in seasons:
            model     = season + str(y)
            if season == 'zomer':
                period  = slice(str(y) + '-04-15',str(y) + '-10-14')
            else:
                period  = slice(str(y) + '-10-15',str(y+1) + '-04-14')
            print(model)
            path      = os.path.join(general.model_dir,'DFM_map_' + season + str(y) + '.nc')
            xds_y     = xr.open_dataset(path).sel(time=period)
            if season + str(y) != 'winter2023': xds_y     = xds_y.sel(time=period)
            xds_yQ    = xds_y[['mesh1d_q1','mesh1d_edge_branch']].groupby('mesh1d_edge_branch').mean().drop_vars('mesh1d_edge_branch')
            xds_yQ    = xds_yQ.transpose("time", "mesh1d_edge_branch")
            xds_y     = xds_yQ[['mesh1d_q1']]
            if model == 'zomer2014':
                xds_var   = xds_y
            else:
                xds_var   = xr.concat([xds_var,xds_y],dim='time')
            xds_y.close()    
            del xds_y
    xds_var.to_netcdf(general.path_Q)
    xds_var.close()
    
   

if __name__ == "__main__":
    
    # Define variables
    fn        = 'DFM_map.nc'  
    
    # prepare data
    if not os.path.exists(general.path_base):
        print('Extract base data')
        get_base_data()
    if not os.path.exists(general.path_WL):
        print('Extract & merge WL data')
        get_WL_data()
    if not os.path.exists(general.path_Q):
        print('Extract & merge Q data')
        get_Q_data()
    if not os.path.exists(general.path_all_WL):
        print('Create merged dataset')
        xds_base = xr.open_dataset(general.path_base)
        xds_WL   = xr.open_dataset(general.path_WL)    
        xds_all  = xr.merge([xds_base,xds_WL])   
        xds_all.to_netcdf(general.path_all_WL)
        xds_WL.close()
        xds_base.close()
        os.remove(general.path_WL)
    if not os.path.exists(general.path_all_Q):
        print('Create merged dataset')
        xds_base = xr.open_dataset(general.path_base)
        xds_Q    = xr.open_dataset(general.path_Q)    
        xds_all  = xr.merge([xds_base,xds_Q])   
        xds_all.to_netcdf(general.path_all_Q)
        xds_base.close()
        xds_Q.close()
        os.remove(general.path_base)
        os.remove(general.path_Q)
        
    