# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 13/04/2026
virtual environment used: geo-env
"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import get_coordinates

from shapely.geometry import Point
from shapely import get_coordinates

import warnings
warnings.filterwarnings("ignore")


class general():
    path_all            = r'D:/workingdir/3_Output/netcdf/DFM_all_Waterstand.nc'
    path_WLstats        = r'D:/workingdir/3_Output/gpkg/Waterstandstatistieken.gpkg'

def get_WLstats(xds):
    var_dat   = 'mesh1d_s1' # WL
    
    # transform to shapefile
    lmax       = len(xds['mesh1d_nNodes'])
    ID         = np.array([i.strip().replace('hdsr_wa_','') for i in xds['mesh1d_node_id'].values.astype(str)])    
    index      = np.arange(0,lmax)
    x_coord    = xds['mesh1d_node_x'].values
    y_coord    = xds['mesh1d_node_y'].values
    df         = pd.DataFrame({'index': index, 'CODE': ID, 'X': x_coord, 'Y': y_coord})
    df['geometry'] = df.apply(lambda x: Point((float(x.X), float(x.Y))), axis=1)
    gdf        = gpd.GeoDataFrame(df, geometry='geometry').set_index('index')
    
    # get statistics
    gdf['WL0']    = xds[var_dat].isel(time=0)
    gdf['WLmax']  = xds[var_dat].max(dim='time').values
    gdf['WLmin']  = xds[var_dat].min(dim='time').values
    gdf['WLgem']  = xds[var_dat].mean(dim='time').values
    gdf['dWLmax'] = gdf['WLmax'] - gdf['WL0']
    gdf['dWLmin'] = gdf['WLmin'] - gdf['WL0']
    gdf['WL_T1']  = xds[var_dat].groupby("time.year").max("time").median("year").values # .isel(time=summer)
    for perc in np.arange(0,1.1,0.1):
        stat =  str(int(perc*100))
        gdf['WL_p' + stat]  = xds[var_dat].quantile(perc,dim='time')
    
    return gdf

def fun_export_shp(gdf):
    # round
    for v in gdf.columns:
        if 'float' in str(gdf[v].dtype):
            gdf[v]  = np.round(gdf[v],2)
        
    # save
    gdf.crs    = 'EPSG:28992'
    gdf.to_file(general.path_WLstats)                



if __name__ == "__main__":
    
    # load data
    print('Loading data...')
    xds = xr.open_dataset(general.path_all)
    
    # get statistics
    print('Calculating statistics...')
    gdf = get_WLstats(xds)
        
    # export shapefile
    fun_export_shp(gdf)
       
    xds.close()
    print('Done!')