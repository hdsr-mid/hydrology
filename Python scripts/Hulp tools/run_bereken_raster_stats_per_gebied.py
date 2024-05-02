# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:27:11 2023

@author: PetraH

Bepaal x-percentiel hoogte
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
import rasterio
from rasterio.mask import mask
from rasterstats import zonal_stats
from shapely.geometry import mapping

# Geef id van de zone op waarbinnen de percentages berekend moet worden
zonalid = 'OBJECTID'

class config():
    # define all paths here    
    root   = os.getcwd()
    
    # Paths: Input files
    afwatering = os.path.join(root, "afwateringseenheden_2020.shp")
    water      = os.path.join(root,'BR_VS_Insteekvlak.shp')
    topo       = os.path.join(root,'AHN3.tif')
    
# functie om percentielen van hoogtemodel te berekenen
def compute_zonal_stats(fp_vector, fp_raster, stats):
    """ Zonal stats voor AHN"""
    gdf = gpd.read_file(fp_vector)

    stats = zonal_stats(
        gdf,
        fp_raster,
        stats=stats
    )

    gdf_sel = gdf[[zonalid, 'geometry']]
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

if __name__ == "__main__":
    shp_afwatering = gpd.read_file(config.afwatering)
    shp_water      = gpd.read_file(config.water)
    topo           = rasterio.open(config.topo)
    stats          = ['percentile_10','median', 'count']
    
    nmax = len(shp_afwatering)
    stat1= np.zeros(nmax) - 999    
    stat2= np.zeros(nmax) - 999    
    for i in tqdm(range(0,nmax)):   
        fn_shp = config.afwatering.replace('.shp','_'+str(i)+'.shp')
        fn_tif = config.topo.replace('.tif','_'+str(i)+'.tif')
        
        row = shp_afwatering.iloc[[i]]
        geometries = row.geometry.apply(mapping)
        
        # clip water
        # selecteer BGT-water in zone peilgebied
        gdf_differ_water = gpd.overlay(row, shp_water[['geometry']], how='difference')
        gdf_differ_water.to_file(fn_shp)
        
        # clip hoogtemodel op geselecteerd gebied
        cropped_ahn, out_transform = mask(topo, geometries, crop=True)
        out_meta = get_updated_meta(topo, cropped_ahn, out_transform)            
        with rasterio.open(fn_tif, 'w', **out_meta) as dest:
            dest.write(cropped_ahn)
            
        gdf_zonal_stats = compute_zonal_stats(fn_shp, fn_tif, stats)
        if gdf_zonal_stats['count'].values != 0:
            stat1[i] = float(gdf_zonal_stats['percentile_10'].values)        
            stat2[i] = float(gdf_zonal_stats['median'].values)        
        
        fn = [os.path.join(os.getcwd(),f) for f in os.listdir(os.getcwd())]
        fn = [f for f in fn if f.startswith(fn_shp.replace('.shp',''))]
        os.remove(fn_tif)
        for f in fn: os.remove(f)
    shp_afwatering['topo_10p'] = stat1
    shp_afwatering['topo_50p'] = stat2
    
    shp_afwatering.to_file(config.afwatering.replace('.shp','_stat.shp'))