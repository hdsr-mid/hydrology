# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 12:13:55 2024

@author: PetraH
"""
''' Add the code and smallest distance to a Hydro object to the profiles'''

import os
import numpy as np
import geopandas as gpd
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")

class paths():
    # Input data
    root          = os.getcwd()    
    shp_HO        = os.path.join(root,'GIS','BR_VS_HydroObject.shp')
    shp_points    = os.path.join(root,'GIS','points.shp') # script based on file with profilepoints
    
    # output data
    shp_out       = os.path.join(root,'GIS','points.gpkg')

if __name__ == "__main__":    
    # Input data
    shp_HO     = gpd.read_file(paths.shp_HO)
    shp_HO     = shp_HO[['CODE','CATEGORIEO','geometry']]
    shp_points = gpd.read_file(paths.shp_points)    
    shp_points.to_file(paths.shp_out, driver="GPKG",layer = 'points_orig')
    
    # Group points based on the column 'profielmet'
    grouped = shp_points.groupby('pm_id')
    
    # Loop through each profielmet-code
    for name, group in tqdm(grouped):    
        # Spatial join
        gdf_joined = gpd.sjoin_nearest(group, shp_HO, distance_col="distances")
        
        # get smallest distance
        distance = np.min(gdf_joined.distances.values)
        
        # add info to gdf        
        shp_points.loc[(shp_points['pm_id'] == name), 'HO_ID']      = gdf_joined.CODE.values[0]
        shp_points.loc[(shp_points['pm_id'] == name), 'HO_cat']     = gdf_joined.CATEGORIEO.values[0]
        shp_points.loc[(shp_points['pm_id'] == name), 'HO_afstd']   = distance
        
    # Save file
    shp_points.to_file(paths.shp_out, driver="GPKG",layer = 'points_incl_HO')
    