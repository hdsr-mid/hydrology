# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 12:15:19 2024

@author: phulsman

In deze script worden eindpunten van watergangen geindentificeerd die niet met elkaar verbonden zijn.

"""

import os
import numpy as np
import geopandas as gpd
from tqdm import tqdm
from shapely.ops import transform
import shapely

class paths():
    root = os.getcwd()
    
    shp_watergangen = os.path.join(root,'GIS_input','br_vs_hydroobject.shp')
    shp_out         = os.path.join(root,'GIS_output','points_niet_verbonden.shp')

if __name__ == "__main__":
        
        # input data
        col_watergangen = ['CODE','CATEGORIEO','WS_MODEL','geometry']
        shp_watergangen = gpd.read_file(paths.shp_watergangen)[col_watergangen]  

        # get endpoints
        gdf             = gpd.GeoDataFrame()
        gdf["CODE"]     = np.append(shp_watergangen.CODE.values,shp_watergangen.CODE.values)
        gdf["geometry"] = np.append(shapely.get_point(shp_watergangen.geometry, 0),shapely.get_point(shp_watergangen.geometry, -1))
        gdf             = gdf.drop_duplicates(subset='geometry', keep=False)
        gdf["dist_p"]   = np.array([gdf.iloc[x].geometry.distance(gdf[gdf.CODE != gdf.iloc[x].CODE].geometry).min() for x in tqdm(range(0,len(gdf)))])        
        gdf["dist_l"]   = np.array([gdf.iloc[x].geometry.distance(shp_watergangen[shp_watergangen.CODE != gdf.iloc[x].CODE].geometry).min() for x in tqdm(range(0,len(gdf)))])
        gdf             = gdf[gdf.dist_p>0.001]        
        gdf             = gdf[(gdf.dist_l<1) | (gdf.dist_p<1)]        
        gdf             = gdf.set_geometry('geometry')
        gdf.crs         = 'EPSG:28992'
        gdf.to_file(paths.shp_out)
        
        # verwijder punten dat direct op een niet-primaire watergang ligt (afstand = 0)
        shp_non_primair = shp_watergangen[shp_watergangen['CATEGORIEO']!=1]
        ind = [x for x in tqdm(range(0,len(gdf))) if gdf.iloc[x].geometry.distance(shp_non_primair.geometry).min() > 0]
        gdf = gdf.iloc[ind]
        gdf.to_file(paths.shp_out.replace('.shp','_sel.shp'))
        
        # verwijder punten dicht bij een lijn (afstand < 0.001 m)
        gdf             = gdf[gdf.dist_l>0.001]        
        gdf.to_file(paths.shp_out.replace('.shp','_final.shp'))
        