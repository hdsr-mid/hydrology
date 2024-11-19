# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 13:15:19 2023

@author: PetraH

In deze script wordt de inundatiediepte berekend o.b.v. een bepaald waterniveau per peilgebied.
Dit waterniveau is gelijk aan het maximum gemodelleerd waterniveau voor een peilgbied.

Als referentie is de 10% percentiel hoogte gebruikt. 
Hier wordt dus de waterdiepte boven het toegestaan waterniveau berekend.


"""

import os
import numpy as np
import pandas as pd
import rasterio
from tqdm import tqdm
import geopandas as gpd
from rasterio.mask import mask
from rasterio.merge import merge
from shapely.geometry import mapping
import warnings
warnings.filterwarnings("ignore")

# id of the zone over which the statistics are calculated
zonalid = 'CODE'

class paths():
    root = os.getcwd()
    # Input files
    tif_AHN         = os.path.join(root, 'Input',"AHN4_ruw_clip_prep.tif") # zie workflow "AHN burn water"
    shp_peil        = os.path.join(root, 'Input','BR_VS_Peilgebieden_WL.shp') # zie workflow "Sobek 2 shapefiles"
    shp_insteek     = os.path.join(root,'Input','BR_VS_Insteekvlak.shp')
    
    # Temp files
    raster_topo       = os.path.join(root, 'Output', "temp_AHN4_clip_peilgebied.tif")
    shp_droognat      = os.path.join(root, 'Output','temp_shp_droognat.shp')   
    
    # Output files
    excel_file       = os.path.join(root, 'Output',"MM_INNUN.txt")
    raster_WD        = os.path.join(root, 'Output',"waterdepth.tif")
    
def fun_inundation(idx, water_level, polygon, inun_type):  
    
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
    raster_AHN  = rasterio.open(paths.tif_AHN)   
    polygons    = gpd.read_file(paths.shp_peil)  
    shp_insteek = gpd.read_file(paths.shp_insteek)
    # shp_insteek = shp_insteek[shp_insteek['CATEGORIEO']==1]
    
    # Calculate polygon areas in square meters
    polygons['area_m2'] = polygons['geometry'].area
        
    # Initialize an empty DataFrame to store the results
    df_volume = pd.DataFrame(columns=['code', 'water_level', 'volume_m3','volume_mm'])
    
    # Iterate over each polygon
    for idx in tqdm(range(0,len(polygons))):
        
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
        fun_inundation(idx, water_level, polygon, inun_type)            
        
        # Estimate volume
        raster_WD  = rasterio.open(paths.raster_WD.replace('.tif',f"_{idx}.tif"))
        data       = raster_WD.read(1)
        pos_data   = np.where(data > 0, data, 0)  # Only consider positive inundation        
        resolution = abs(raster_WD.transform[0] * raster_WD.transform[4])  # Calculate pixel area
        volume     = np.sum(pos_data/100) * resolution
        raster_WD.close()
        
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
    
    # Use the to_excel method to export the DataFrame to Excel
    df_volume.to_csv(paths.excel_file, index=False)  # Set index to False if you don't want to export the row numbers as a column
    

    