# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:57:52 2024

@author: PetraH

Plot verschillende onderdelen van DHydro:
    - .pliz         fixed weir (hoogtelijnen)
    - .bc           boundary condition
    - _his.nc       output time series 
    - _map.nc       output gridded
    
virtual environment: dsd_env

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import contextily as ctx
import xugrid as xu

if __name__ == "__main__":
    
    ''' .pliz file (fixed weir/hoogtelijnen) ''' 
    fn = r'DATA/objects.pli_fxw.pliz'
    polyfile_object = hcdfm.PolyFile(fn)
    gdf_polyfile = dfmt.PolyFile_to_geodataframe_linestrings(polyfile_object,crs='EPSG:28992')    
    ax = gdf_polyfile.plot()
    ax.set_aspect('equal')
    plt.savefig(r'figs_output/fig_pliz.png')
    plt.close()
    
    ''' .bc file (boundary condition) ''' 
    fn  = r'DATA/rhine.bc'
    var = 'lateral_discharge'
    forcingmodel_object = hcdfm.ForcingModel(fn)
    fig, ax = plt.subplots()
    for iFO, forcingobj in enumerate(forcingmodel_object.forcing):
        forcing_xr = dfmt.forcinglike_to_Dataset(forcingobj, convertnan=True)
        forcing_xr[var].plot(ax=ax, label=forcing_xr[var].attrs['locationname'], linewidth=0.8)
    ax.legend(loc=1)    
    plt.savefig(r'figs_output/fig_bc.png')
    plt.close()
    
    