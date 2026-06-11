# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Clip rr/.nc op basis van id_keep
    

"""
import os
import xarray as xr

def clip_rr_netcdf(path, id_keep):
    # input data
    files= [f for f in os.listdir(path) if f.endswith('.nc')]
    
    # iterate through files
    for f in files:
        # load data
        fn_in  = f
        fn_out = fn_in.replace('.nc','_sel.nc')
        xds    = xr.open_dataset(os.path.join(path,fn_in))
        print('clip ' + fn_in)
        
        # get location_names
        id_f    = [i.strip() for i in xds.locations.values.astype(str)]
        i_sel   = [i for i in range(len(id_f)) if id_f[i] in id_keep]
        
        # select
        xds_sel = xds.isel(id=i_sel)
        
        # save file
        xds_sel.to_netcdf(os.path.join(path,fn_out))
        
        # remove old file, rename new
        xds.close(), xds_sel.close()
        del xds, xds_sel
        os.remove(os.path.join(path,fn_in))
        os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
        
    
