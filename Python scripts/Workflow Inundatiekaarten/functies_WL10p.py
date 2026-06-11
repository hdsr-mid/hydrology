# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:32:57 2023

@author: PetraH
"""

import os
import numpy as np
import pandas as pd
import rasterio
from tqdm import tqdm
import geopandas as gpd
from rasterio.mask import mask
import config
from shapely.geometry import mapping
from rasterstats import zonal_stats

# Geef id van de zone op waarbinnen de percentages berekend moet worden
CODE = config.CODE


# functie om percentielen van hoogtemodel te berekenen
def compute_zonal_stats(fp_vector, fp_raster, stats):
    """ Zonal stats voor AHN"""
    gdf = gpd.read_file(fp_vector)

    stats = zonal_stats(
        gdf,
        fp_raster,
        stats=stats
    )
    gdf_sel = gdf[[CODE, 'geometry']]
    df_concat = pd.concat((gdf_sel, pd.DataFrame(stats)), axis=1)
    return df_concat

#    functie om uitsnede bodemhoogtemodel op te slaan
def get_updated_meta(f, out_img, out_transform):
    """ Raster meta-informatie"""
    out_meta = f.meta.copy()
    out_meta['transform'] = out_transform
    out_meta['compress'] = 'deflate'
    out_meta['height'] = out_img.shape[1]
    out_meta['width'] = out_img.shape[2]
    #out_meta['crs'] = {'init': 'epsg:28992'}
    out_meta['crs'] = None
    return out_meta
    

def polygon_add_10p(shp_poly):
    # paths
    raster_topo          = config.raster_prep
    shp_insteek          = config.shp_insteek
    
    # Input data
    shp_water      = gpd.read_file(shp_insteek)
    topo           = rasterio.open(raster_topo)
    stats          = ['percentile_10','median', 'count']
        
    # empty array
    MV_10p = np.zeros(len(shp_poly)) - 999
    
    # Iterate over each polygon
    for i in tqdm(range(0,len(shp_poly)),'Estimate MV 10p'):
        # Input
        fn_shp         = config.shp_poly_10p.replace('.shp','_'+str(i)+'.shp')
        fn_tif         = raster_topo.replace('.tif','_'+str(i)+'.tif')    
        row            = shp_poly.iloc[[i]]
        geometries     = row.geometry.apply(mapping)
        
        # clip water
        # selecteer BGT-water in zone peilgebied
        gdf_differ_water = gpd.overlay(row, shp_water[['geometry']], how='difference')
        gdf_differ_water.to_file(fn_shp)
        
        # clip hoogtemodel op geselecteerde peilgebied
        cropped_ahn, out_transform = mask(topo, geometries, crop=True)
        out_meta = get_updated_meta(topo, cropped_ahn, out_transform)            
        with rasterio.open(fn_tif, 'w', **out_meta) as dest:
            dest.write(cropped_ahn)
            
        gdf_zonal_stats = compute_zonal_stats(fn_shp, fn_tif, stats)
        if gdf_zonal_stats['count'].values != 0:
            stat      = float(gdf_zonal_stats['percentile_10'].values)
            MV_10p[i] = stat
            
        fn = [os.path.join(os.getcwd(),f) for f in os.listdir(os.getcwd()) if f.startswith(fn_shp.replace('.shp',''))]
        os.remove(fn_tif)
        for f in fn: os.remove(f)
        
    shp_poly['MV_10p'] = MV_10p
    
    return shp_poly
            
def fun(shp_poly):
    shp_poly_10p = shp_poly.copy()
    
    # add 10%-percentile to polygon
    shp_poly_10p = polygon_add_10p(shp_poly_10p)
    
    return shp_poly_10p