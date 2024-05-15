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

class paths():
    root = os.getcwd()
    # Input files
    raster_prep     = os.path.join(root, 'Output',"AHN3_ruw_clip_prep.tif")
    peil_clip       = os.path.join(root, 'Output','BR-VS_Peilgebieden_clip_WL_AHN10p.shp')
    shp_insteek     = os.path.join(root,'Input','BR_VS_Insteekvlak.shp')
    
    # Temp files
    raster_topo       = os.path.join(root, 'Output', "AHN3_clip_peilgebied.tif")
    raster_WD_binary  = os.path.join(root, 'Output', "waterdepth_binary.tif")
    shp_droognat      = os.path.join(root, 'Output','shp_droognat.shp')   
    
    # Output files
    excel_file       = os.path.join(root, 'Output',"MM_INNUN.txt")
    raster_WD        = os.path.join(root, 'Output',"waterdepth.tif")
    raster_inun      = os.path.join(root, 'Output',"inundationdepth.tif")
    shp_out_inun     = os.path.join(root, 'Output','inundatie.shp') 
    
def fun_inundation_p1(idx, water_level, polygon, inun_type):  
    ''' Step 1: area below waterlevel -> waterdepth'''
    # Open height data
    raster     = rasterio.open(paths.raster_topo)
    meta       = raster.meta.copy()
    meta.update(compress='lzw')
    height_data = raster.read(1)
        
    if inun_type != 'total': 
        # Estimate 10%-percentile height within polygon
        height_10p  = polygon['AHN_10p'].values
        
        # Edit height data (remove area below 10%-percentile)        
        height_filled   = np.where(height_data < height_10p, height_10p, height_data)
        height_filled   = np.where(height_data < -900, height_data, height_filled) # no-data pixels need to remain no-data pixels
        height_data     = height_filled
    
    # Add water and calculate waterdepth
    waterdepth = water_level - height_data
    waterdepth = np.where(np.logical_or(waterdepth < 0, height_data < -500), 0, waterdepth)
    waterdepth = np.round(100*waterdepth,0) # [m] -> [cm]        
        
    # Create a new raster file with waterdepth values
    with rasterio.open(
        paths.raster_WD.replace('.tif',f"_{idx}.tif"),
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
    
    
def fun_inundation_p2(polygon, idx, shp_insteek):
    ''' Step 2: area below waterlevel which is connected to a waterbody -> inundationdepth'''
    # classify inundation map (0: dry, 1: wet)
    raster   = gdal.Open(paths.raster_WD.replace('.tif',f"_{idx}.tif"))
    band     = raster.GetRasterBand(1)
    arr      = band.ReadAsArray()
    arr      = np.where(arr> 0, 1, 0)
    gdal_array.SaveArray(arr.astype("int16"), paths.raster_WD_binary, "GTIFF", raster)
    del raster
    
    # polygonize
    raster   = gdal.Open(paths.raster_WD_binary)
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
    shp        = gpd.read_file(paths.shp_out_inun)
    raster     = rasterio.open(paths.raster_WD.replace('.tif',f"_{idx}.tif"))
    out_meta   = raster.profile
    geom_value = ((geom,value) for geom, value in zip(shp.geometry, shp['Inundatie']))
    rasterized = features.rasterize(geom_value,
                        out_shape = raster.shape,
                        transform = raster.transform,
                        all_touched = False,
                        fill = 0,   # background value
                        merge_alg = MergeAlg.replace,
                        dtype = "int16")
    data = np.where(rasterized == 1, raster.read(1), 0)
    with rasterio.open(paths.raster_inun.replace('.tif',f"_{idx}.tif"), 'w', **out_meta) as dest:
        dest.write(data,indexes=1)
    raster.close()
    
    # remove file
    del shp, shp_clip
    os.remove(paths.shp_droognat)
    os.remove(paths.raster_WD_binary)

def fun_merge_tif(tif_file, polygons):    
    tif_files = [tif_file.replace('.tif',f"_{idx}.tif") for idx, polygon in polygons.iterrows()]
    inundation_map, out_trans = merge(tif_files)
    src      = rasterio.open(tif_files[0])
    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff","height": inundation_map.shape[1],"width": inundation_map.shape[2],"transform": out_trans,})
    out_meta.update({"nodata": 0})
    with rasterio.open(tif_file, "w", **out_meta) as dest:
        dest.write(inundation_map)
    src.close()
    for f in tif_files: os.remove(f)    

def fun_merge_shp(shp_file, polygons):
    path = [shp_file.replace('.shp',f"_{idx}.shp") for idx, polygon in polygons.iterrows()]
    shp  = gpd.GeoDataFrame(pd.concat([gpd.read_file(i) for i in path],ignore_index=True), crs=gpd.read_file(path[0]).crs)
    shp.to_file(paths.shp_file.replace('.shp','.gpkg'), driver="GPKG",layer = 'inundatie')
       
    
if __name__ == "__main__":
    # Load data
    raster_AHN  = rasterio.open(paths.raster_prep)   
    polygons    = gpd.read_file(paths.peil_clip)  
    shp_insteek = gpd.read_file(paths.shp_insteek)
    # shp_insteek = shp_insteek[shp_insteek['CATEGORIEO']==1]
    
    # Calculate polygon areas in square meters
    polygons['area_m2'] = polygons['geometry'].area
        
    # Initialize an empty DataFrame to store the results
    df_volume = pd.DataFrame(columns=['code', 'water_level', 'volume_m3','volume_mm'])
    
    # Iterate over each polygon
    for idx in range(0,len(polygons)):
        print(idx,' of ',len(polygons))
        # Extract the geometry of the polygon
        polygon      = polygons[idx:idx+1]
        geometry     = polygon.geometry.apply(mapping)
        
        # Extract water level
        water_level = float(polygon['WLmax'].values)
        
        # Clip the raster with the polygon geometry
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
        
        # Estimate inundation
        inun_type   = 'total' # 'not_total'
        fun_inundation_p1(idx, water_level, polygon, inun_type)            
        fun_inundation_p2(polygon, idx, shp_insteek)            
        
        # Estimate volume
        raster_inun= rasterio.open(paths.raster_inun.replace('.tif',f"_{idx}.tif"))
        data       = raster_inun.read(1)
        pos_data   = np.where(data > 0, data, 0)  # Only consider positive inundation        
        resolution = abs(raster_inun.transform[0] * raster_inun.transform[4])  # Calculate pixel area
        volume     = np.sum(pos_data/100) * resolution
        raster_inun.close()
        
        if volume is not None:
            # Calculate mm only if volume_step is not None
            volume_mm = float(volume / polygon['area_m2'].values * 1000)
            
            # Append to the DataFrame
            df_volume = pd.concat([df_volume, pd.DataFrame([{
                'code': polygons['CODE'][idx],
                'water_level': water_level,
                'volume_m3': np.round(volume,2),
                'volume_mm': np.round(volume_mm,0),                    
            }])], ignore_index=True)
        else:
            print(f"Skipping water level {water_level} for {polygons['CODE'][idx]} due to an error in processing.")
            
    
    # Merge tif files
    fun_merge_tif(paths.raster_WD, polygons)
    fun_merge_tif(paths.raster_inun, polygons)
    
    # Use the to_excel method to export the DataFrame to Excel
    df_volume.to_csv(paths.excel_file, index=False)  # Set index to False if you don't want to export the row numbers as a column
    

    