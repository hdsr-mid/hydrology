# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 13/04/2026
virtual environment used: geo-env
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from shapely.geometry import MultiLineString
from shapely import get_coordinates
import warnings
warnings.filterwarnings("ignore")

class general():
    path_all            = r'D:/workingdir/3_Output/netcdf/DFM_all_Debiet.nc'
    shp_watergangen     = r'D:/workingdir/1_InputData/Branches.gpkg'
    path_Qstats         = r'D:/workingdir/3_Output/gpkg/Debietstatistieken.gpkg'
    
def getLineCoords(row, geom, coord_type):
    if isinstance(row[geom], MultiLineString):
        if coord_type == 'x':
            return get_coordinates(row[geom])[:,0]
        elif coord_type == 'y':
            return get_coordinates(row[geom])[:,1]
    else:
        if coord_type == 'x':
            return list( row[geom].coords.xy[0] )
        elif coord_type == 'y':
            return list( row[geom].coords.xy[1] )
        
def fun_Qstats(xds):
    gdf = gpd.read_file(general.shp_watergangen)
    gdf['CODE'] = np.array([i.strip().replace('hdsr_wa_','') for i in gdf['Name'].values])
    gdf = gdf.drop(["Name"], axis=1)
    
    ID  = np.array([i.decode().strip().replace('hdsr_wa_','') for i in xds['network_branch_id'].values])
    idx1 = np.array([i for i in range(0,len(gdf)) if gdf.CODE.values[i] in ID])
    idx2 = np.array([i for i in range(0,len(gdf)) if gdf.CODE.values[i] not in ID])
    gdf1 = gdf.iloc[idx1]
    gdf2 = gdf.iloc[idx2]
    order = np.array([int(np.where(g==ID)[0]) for g in gdf1.CODE.values])
    
    # get sumer days
    months = [pd.to_datetime(t).month for t in xds.time.values]
    summer = np.array([i for i in range(0,len(months)) if months[i] in [5,6,7,8,9]])
    
    # get statistics
    var_dat       = 'mesh1d_q1' # Q
    gdf1['Q_gem']  = xds[var_dat].mean(dim='time').values[order]
    gdf1['Q_T1']   = xds[var_dat].groupby("time.year").max("time").median("year").values[order] # .isel(time=summer)
    gdf2['Q_gem']  = -999
    gdf2['Q_T1']   = -999
    for perc in np.arange(0,1.1,0.1):
        stat =  str(int(perc*100))
        gdf1['Q_p' + stat]  = xds[var_dat].quantile(perc,dim='time').values[order]
        gdf2['Q_p' + stat]  = -999
    gdf = pd.concat([gdf1,gdf2])
    
    gdf = pd.concat([gdf1,gdf2])
    return gdf

def fun_export_shp(gdf):
    # round
    for v in gdf.columns:
        if 'float' in str(gdf[v].dtype):
            gdf[v]  = np.round(gdf[v],2)
        
    # save
    gdf.to_file(general.path_Qstats)



if __name__ == "__main__":
    
    # load data
    print('Loading data...')
    xds_all = xr.open_dataset(general.path_all)
    
    # get statistics
    print('Calculating statistics...')
    gdf = fun_Qstats(xds_all)
        
    # export shapefile
    fun_export_shp(gdf)
    
    xds_all.close()
    print('Done!')