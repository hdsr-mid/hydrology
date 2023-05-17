# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:04:32 2023

@author: PetraH
"""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

def fun_plot(xds):
    # select data
    data = xds['stations']
    
    # select time step by index    
    data = data.isel(time=0)
    
    # transpose data
    data = data.transpose('y','x')
    
    # plot
    plt.figure(figsize=(10,6))
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplot2grid((1,1), (0,0))
    
    plt.imshow(data,cmap='seismic')
    plt.colorbar(label='data')
    
    plt.savefig('fig.png',dpi=300)
    
if __name__ == '__main__':
    # file name
    fn = 'rain.nc'
    
    # Read netcdf file
    xds = xr.open_dataset(fn)
    print(xds)
    
    # loop over variables
    for var in list(xds.head()):
        print(var)
        
        # select variable
        xds_var = xds[[var]]
        
        # save netcdf file
        encoding = {var: dict(zlib=True, complevel=9) for var in xds_var.data_vars}
        xds_var.to_netcdf(fn.replace('.nc','_' + var+'.nc'),encoding=encoding)
    
        
    # Plot data for specific variable and time step
    fun_plot(xds)
    
    
    