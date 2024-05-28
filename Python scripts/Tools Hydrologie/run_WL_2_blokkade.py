# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 13:15:19 2023

@author: PetraH
"""

import os
import numpy as np
import pandas as pd
import rasterio
from tqdm import tqdm
import geopandas as gpd
from rasterio.mask import mask
from rasterio.merge import merge
from osgeo import gdal, ogr, gdal_array
from rasterio import features
from rasterio.enums import MergeAlg
from shapely.geometry import mapping
import warnings
warnings.filterwarnings("ignore")

'''The inundatie is estimated in two steps:
    1. area below waterlevel -> waterdepth
    2. area below waterlevel which is connected to a waterbody -> inundationdepth '''

# id of the zone over which the statistics are calculated
zonalid = 'CODE'
threshold = 20 # a waterdepth of more than this value in [cm] is considered a blokkage


class paths():
    root = os.getcwd()
    
    # Input files
    raster_prep     = os.path.join(root, 'Output',"AHN3_ruw_clip_prep.tif")
    peil_clip       = os.path.join(root, 'Output','BR-VS_Peilgebieden_clip_edited_10p_60u_LWS_GHG_Wlmax.shp')
    shp_insteek     = os.path.join(root,'Input','BR_VS_Insteekvlak.shp')
    
    # Temp files
    raster_topo      = os.path.join(root, 'Output','temp', "AHN3_clip_peilgebied.tif")
    raster_binary    = os.path.join(root, 'Output','temp', "waterdepth_binary.tif")
    shp_droognat     = os.path.join(root, 'Output','temp','shp_droognat.shp')   
    raster_WD        = os.path.join(root, 'Output','temp',"waterdepth.tif")
    raster_inun      = os.path.join(root, 'Output','temp',"inundationdepth.tif")
    shp_out_inun     = os.path.join(root, 'Output','temp','inundatie.shp') 
    shp_blok         = os.path.join(root, 'Output','temp','water_blokkade.shp') 
    
    # Output files
    shp_out          = os.path.join(root, 'Output','water_blokkade.shp') 

def fun_remove_shp(shp_path):
    dir_name    = os.path.dirname(shp_path)
    file_prefix = os.path.basename(shp_path).replace('.shp','')
    files       = [os.path.join(dir_name,f) for f in os.listdir(dir_name) if file_prefix in f]
    
    for f in files: os.remove(f)
    
    
def fun_clip_AHN(raster_AHN,polygon):
    # Get geometry
    geometry     = polygon.geometry.apply(mapping)
    
    # Clip
    topo_clip, out_transform = mask(raster_AHN, geometry, crop=True)
    
    # Update the metadata
    topo_meta = raster_AHN.meta.copy()
    topo_meta.update({
        "driver": "GTiff",
        "height": topo_clip.shape[1],
        "width": topo_clip.shape[2],
        "transform": out_transform,
    })
    
    # Write the clipped raster to a new file
    with rasterio.open(paths.raster_topo, "w", **topo_meta) as dest:
        dest.write(topo_clip)
    
def fun_inundation_p1(water_level, polygon):  
    ''' Step 1: area below waterlevel -> waterdepth'''
    # Open height data
    raster     = rasterio.open(paths.raster_topo)
    meta       = raster.meta.copy()
    meta.update(compress='lzw')
    height_data = raster.read(1)
        
    # Add water and calculate waterdepth
    waterdepth = water_level - height_data
    waterdepth = np.where(np.logical_or(waterdepth < 0, height_data < -500), 0, waterdepth)
    waterdepth = np.round(100*waterdepth,0) # [m] -> [cm]        
        
    # Create a new raster file with waterdepth values
    with rasterio.open(
        paths.raster_WD,
        'w',
        driver='GTiff',
        height=raster.height,
        width=raster.width,
        count=1,
        compress='lzw',
        dtype='int16',
        crs=raster.crs,
        transform=raster.transform,
        nodata = 0,
    ) as output_dataset:
        output_dataset.write(waterdepth, 1)
     
    raster.close()
    del raster
     
    # Remove files
    os.remove(paths.raster_topo)
    
def fun_inundation_p2(polygon, shp_insteek):
    ''' Step 2: area below waterlevel which is connected to a waterbody -> inundationdepth'''
    # classify inundation map (0: dry, 1: wet)
    raster   = gdal.Open(paths.raster_WD)
    band     = raster.GetRasterBand(1)
    arr      = band.ReadAsArray()
    arr      = np.where(arr> 0, 1, 0)
    gdal_array.SaveArray(arr.astype("int16"), paths.raster_binary, "GTIFF", raster)
    del raster
    
    # polygonize
    raster   = gdal.Open(paths.raster_binary)
    band     = raster.GetRasterBand(1)
    drv      = ogr.GetDriverByName('ESRI Shapefile')
    outfile  = drv.CreateDataSource(paths.shp_droognat)
    outlayer = outfile.CreateLayer('polygonized raster', srs = None )
    newField = ogr.FieldDefn('Inundatie', ogr.OFTReal)
    outlayer.CreateField(newField)
    gdal.Polygonize(band, None, outlayer, 0, [])
    outfile  = None
    del raster
    
    # clip
    shp      = gpd.read_file(paths.shp_droognat)
    shp.crs  = polygon.crs 
    shp_clip = gpd.clip(shp, polygon)
    shp_clip.crs='EPSG:28992'
    shp_clip['area'] = shp_clip.area        
    
    # assign value based on shp_insteek
    shp       = shp_clip.copy()
    shp       = shp.loc[(shp['Inundatie']!=0) & (shp['Inundatie']!=-999)].reset_index()       # droge gebieden skippen      
    intersect = np.zeros(len(shp))    
    for index, orig in tqdm(shp.iterrows(),str(len(shp))):               
        count = orig['geometry'].intersects(shp_insteek['geometry']).sum()
        if count > 0:
            intersect[index] = 1
    shp['intersect'] = intersect
    shp = shp.loc[shp['intersect']==1]
    shp.to_file(paths.shp_out_inun)
    
    # Rasterize
    raster     = rasterio.open(paths.raster_WD)
    out_meta   = raster.profile
    shp        = gpd.read_file(paths.shp_out_inun)
    if len(shp)>0:
        geom_value = ((geom,value) for geom, value in zip(shp.geometry, shp['Inundatie']))
        rasterized = features.rasterize(geom_value,
                            out_shape = raster.shape,
                            transform = raster.transform,
                            all_touched = False,
                            fill = 0,   # background value
                            merge_alg = MergeAlg.replace,
                            dtype = "int16")
        data = np.where(rasterized == 1, raster.read(1), 0)
    else:
        data = raster.read(1) * 0
        
    with rasterio.open(paths.raster_inun, 'w', **out_meta) as dest:
        dest.write(data,indexes=1)
    raster.close()
    del raster
    
    # remove file
    del shp, shp_clip
    fun_remove_shp(paths.shp_droognat)
    os.remove(paths.raster_binary)
    fun_remove_shp(paths.shp_out_inun)
    os.remove(paths.raster_WD)
    
def fun_blokkade(idx):
    # Identify blokkage
    raster   = gdal.Open(paths.raster_inun)
    band     = raster.GetRasterBand(1)
    arr      = band.ReadAsArray()
    arr      = np.where(arr> threshold, 1, 0)
    gdal_array.SaveArray(arr.astype("int16"), paths.raster_binary, "GTIFF", raster)
    del raster
    
    # polygonize
    raster   = gdal.Open(paths.raster_binary)
    band     = raster.GetRasterBand(1)
    drv      = ogr.GetDriverByName('ESRI Shapefile')
    outfile  = drv.CreateDataSource(paths.shp_blok.replace('.shp',f"_{idx}.shp"))
    outlayer = outfile.CreateLayer('polygonized raster', srs = None )
    newField = ogr.FieldDefn('Inundatie', ogr.OFTReal)
    outlayer.CreateField(newField)
    gdal.Polygonize(band, None, outlayer, 0, [])
    outfile  = None
    del raster
    
    # remove file
    os.remove(paths.raster_inun)
    os.remove(paths.raster_binary)
    
def fun_merge_shp(shp_file, polygons):
    # Merge files
    path = [shp_file.replace('.shp',f"_{idx}.shp") for idx, polygon in polygons.iterrows()]
    shp  = gpd.GeoDataFrame(pd.concat([gpd.read_file(i) for i in path],ignore_index=True), crs=gpd.read_file(path[0]).crs)
    shp  = shp[shp['Inundatie']==1]
    shp.to_file(paths.shp_out)
    
    # Remove files
    for idx in range(0,len(polygons)): 
        fun_remove_shp(shp_file.replace('.shp',f"_{idx}.shp"))
    
if __name__ == "__main__":
    # Load data
    raster_AHN  = rasterio.open(paths.raster_prep)   
    polygons    = gpd.read_file(paths.peil_clip) 
    shp_insteek = gpd.read_file(paths.shp_insteek)
    
    # Calculate polygon areas in square meters
    polygons['area_m2'] = polygons['geometry'].area
        
    # Initialize an empty DataFrame to store the results
    df_volume = pd.DataFrame(columns=['code', 'water_level', 'volume_m3','volume_mm'])
    
    # Iterate over each polygon
    for idx in range(0,len(polygons)):
        print(idx,' of ',len(polygons))
        # Extract the geometry of the polygon
        polygon      = polygons[idx:idx+1]
        
        # Extract water level
        water_level = float(polygon['WLmax'].values)
        
        # Clip the raster with the polygon geometry
        fun_clip_AHN(raster_AHN,polygon)
        
        # Estimate inundation
        fun_inundation_p1(water_level, polygon)            
        fun_inundation_p2(polygon, shp_insteek)            
        
        # Define blokkage
        fun_blokkade(idx)
        
    # Merge files
    fun_merge_shp(paths.shp_blok, polygons)
    
    