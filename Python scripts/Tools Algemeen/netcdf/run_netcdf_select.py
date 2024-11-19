# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:47:08 2024

@author: PetraH
"""

import os
import xarray as xr
import pandas as pd
import datetime


root   = os.getcwd()


if __name__ == "__main__":

    band         = 'Mesh2d_waterdepth'

    
    # Get prepare data
    xds        = xr.open_dataset('DFM_map.nc')
    tinit      = 96 # spin-up time in [hours]
    t_init     = pd.to_datetime(xds.time.values[0]) + datetime.timedelta(hours=tinit)
    xds        = xds.sel(time=slice(t_init, xds.time.values[-1]))
    xds        = xds[[band,'Mesh2d_ucmag','Mesh2d_face_x_bnd','Mesh2d_face_y_bnd']]
    xds.to_netcdf('data.nc')