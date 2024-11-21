# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:57:52 2024

@author: PetraH

Plot verschillende onderdelen van netwerk.nc
    
virtual environment: dsd_env

"""

import xarray as xr
import matplotlib.pyplot as plt
import xugrid as xu

if __name__ == "__main__":
    
    fn = r'DATA/network.nc'
    ds_map = xr.open_dataset(fn)
    ds_map = xu.core.wrap.UgridDataset(ds_map)
    
    variables = list(ds_map.variables.keys())
    
    for var in variables:
        try:
            plt.figure()
            ds_map[var].ugrid.plot()
            plt.savefig(r'figs_netwerk/' + var + '.png')
            plt.close()
        except:
            a=0
        