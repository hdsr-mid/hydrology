# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:28:08 2024

@author: PetraH
"""

import os
import rasterio
import numpy as np
import xarray as xr
import itertools
import geopandas as gpd
import pandas as pd
from osgeo import gdal, ogr, osr, gdal_array
from rasterio import features
from shapely.geometry import shape
from osgeo_utils.gdal_polygonize import gdal_polygonize
from tqdm import tqdm
from rasterio.mask import mask
from rasterio import features
from rasterio.enums import MergeAlg
from shapely.geometry import mapping
import warnings
warnings.filterwarnings("ignore")
import config

CODE = config.CODE

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
    
def main(tif_inundatie):
    
    # Input data
    shp_insteek = gpd.read_file(config.shp_insteek)
    shp_insteek = shp_insteek[shp_insteek['CATEGORIEO']==1]
    shp_peil    = gpd.read_file(config.shp_peil)
    tif_verbonden = tif_inundatie.replace('.tif','_verbonden.tif')
    shp_peil = shp_peil.sort_values(CODE)
    
    # Peil
    for afv in tqdm(range(0,len(shp_peil))):
        # clip by afvoergebied       
        shp_extent   = shp_peil[afv:afv+1]
        raster       = rasterio.open(tif_inundatie)
        geometries   = shp_extent.geometry.apply(mapping)
        
        raster_clip, out_transform = mask(raster, geometries, crop=True)
        out_meta = get_updated_meta(raster, raster_clip, out_transform)            
        with rasterio.open(config.tif_clip, 'w', **out_meta) as dest:
            dest.write(raster_clip)
        raster.close()
        del raster, raster_clip
        
        # classify inundation map (0: dry, 1: wet)
        raster   = gdal.Open(config.tif_clip)
        band     = raster.GetRasterBand(1)
        arr      = band.ReadAsArray()
        arr      = np.where(arr> 0, 1, 0)
        gdal_array.SaveArray(arr.astype("int16"), config.tif_clip_binar, "GTIFF", raster)
        del raster
        
        # polygonize
        raster   = gdal.Open(config.tif_clip_binar)
        band     = raster.GetRasterBand(1)
        drv      = ogr.GetDriverByName('ESRI Shapefile')
        outfile  = drv.CreateDataSource(config.shp_droognat)
        outlayer = outfile.CreateLayer('polygonized raster', srs = None )
        newField = ogr.FieldDefn('Inundatie', ogr.OFTReal)
        outlayer.CreateField(newField)
        gdal.Polygonize(band, None, outlayer, 0, [])
        outfile  = None
        del raster
        
        # clip
        shp      = gpd.read_file(config.shp_droognat)
        shp.crs  = shp_extent.crs 
        shp_clip = gpd.clip(shp, shp_extent)
        shp_clip.crs='EPSG:28992'
        shp_clip['area'] = shp_clip.area        
        
        # assign value based on shp_insteek
        shp       = shp_clip.copy()
        shp       = shp.loc[(shp['Inundatie']!=0) & (shp['Inundatie']!=-999)].reset_index()       # droge gebieden skippen      
        intersect = np.zeros(len(shp))    
        # for index, orig in tqdm(shp.iterrows(),str(len(shp))):               
        for index, orig in shp.iterrows():               
            count = orig['geometry'].intersects(shp_insteek['geometry']).sum()
            if count > 0:
                intersect[index] = 1
        shp['intersect'] = intersect
        shp = shp.loc[shp['intersect']==1]
        shp.to_file(config.shp_out_inun.replace('.shp','_afv' + str(afv) + '.shp'))
        
        # remove file
        del shp, shp_clip
        os.remove(config.shp_droognat)
        os.remove(config.tif_clip)
    
    # merge files
    path = [config.shp_out_inun.replace('.shp','_afv' + str(afv) + '.shp') for afv in range(0,len(shp_peil))]
    shp  = gpd.GeoDataFrame(pd.concat([gpd.read_file(i) for i in path],ignore_index=True), crs=gpd.read_file(path[0]).crs)
    shp.to_file(config.gpk_out, driver="GPKG",layer = 'inundatie')
   
    # select only regions connected to waterway
    # set raster to -999 if shp['Inundatie'] != 1
    raster       = rasterio.open(tif_inundatie)
    shp          = gpd.read_file(config.gpk_out, layer = 'inundatie')  
    out_meta     = raster.profile
    geom_value = ((geom,value) for geom, value in zip(shp.geometry, shp['Inundatie']))
    rasterized = features.rasterize(geom_value,
                        out_shape = raster.shape,
                        transform = raster.transform,
                        all_touched = False,
                        fill = -999,   # background value
                        merge_alg = MergeAlg.replace,
                        dtype = "int16")
    data = np.where(rasterized == 1, raster.read(1), -999)
    with rasterio.open(tif_verbonden, 'w', **out_meta) as dest:
        dest.write(data,indexes=1)
    raster.close()
        