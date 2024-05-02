# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd
import xarray as xr

if __name__ == "__main__":
    temperature = 15 + 8 * np.random.randn(2, 2, 3)
    lon = [[-99.83, -99.32], [-99.79, -99.23]]
    lat = [[42.25, 42.21], [42.63, 42.59]]
    x   = [-99.83, -99.32]
    y   = [-42.25,42.21]
    time = pd.date_range("2014-09-06", periods=3)
    reference_time = pd.Timestamp("2014-09-05")
    
    
    ds0 = xr.Dataset(data_vars={"temperature1":(["x","y","time"],temperature), 
                               "temperature2":(["x","y","time"],temperature)}, 
                    coords={'x':x,'y':y,'time':time}, )
    
    
    
    ds1 = xr.Dataset(data_vars={"temperature1":(["x","y","time"],temperature), 
                               "temperature2":(["x","y","time"],temperature)}, 
                    coords=dict(lon=(["x", "y"], lon),lat=(["x", "y"], lat),time=time), )
    
    
    
    
    ds2 = xr.DataArray(temperature, coords=dict(lon=(["x", "y"], lon),lat=(["x", "y"], lat),time=time), 
                      dims=["x","y","time"],
                      attrs=dict(description="Ambient temperature.",units="degC",),)
    
    print('ds0')                   
    print(ds0)
    print('')
    print('ds1')                   
    print(ds1)
    print('')
    print('ds2')                   
    print(ds2)       
    
    
    encoding = {var: dict(zlib=True, complevel=9) for var in ds0.data_vars}  
    ds0.to_netcdf('file.nc',encoding=encoding)


       