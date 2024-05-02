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
import config

def add_water_to_raster(input_file, output_file, water_level):       
    
    with rasterio.open(input_file) as dataset:
        # Read the raster data as a numpy array
        height_data = dataset.read(1)

        # Add water and calculate inundation
        inundation = water_level - height_data
        inundation = np.where(np.logical_or(inundation < 0, height_data < -500), 0, inundation)
        resolution = abs(dataset.transform[0] * dataset.transform[4])  # Calculate pixel area
        positive_inundation = np.where(inundation > 0, inundation, 0)  # Only consider positive inundation
        volume_step = np.sum(positive_inundation) * resolution
        
        # Create a new raster file with inundation values
        with rasterio.open(
            output_file,
            'w',
            driver='GTiff',
            height=dataset.height,
            width=dataset.width,
            count=1,
            compress='lzw',
            dtype=inundation.dtype,
            crs=dataset.crs,
            transform=dataset.transform,
            nodata = 0,
        ) as output_dataset:
            output_dataset.write(inundation, 1)

    return volume_step

def remove_insteek(extent, shp_insteek, input_file, output_file, insteek_file, shp_insteek_temp):    
   
    # Get shp_insteek
    shp_insteek = gpd.read_file(shp_insteek)
    extent      = gpd.GeoDataFrame([extent])    
    extent.crs  = shp_insteek.crs 
    shp_insteek = gpd.overlay(extent, shp_insteek, how='intersection',keep_geom_type=True)
    shp_insteek['Value'] = 1
    shp_insteek.to_file(shp_insteek_temp)  
    
    # Rasterize shp_insteek
    raster     = rasterio.open(input_file)
    meta       = raster.meta.copy()
    meta.update(compress='lzw')
    with rasterio.open(insteek_file, 'w+', **meta) as out:
        out_arr = out.read(1)
        shapes = ((geom,value) for geom, value in zip(shp_insteek.geometry, shp_insteek.Value))
        burned = rasterio.features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)
    
    # Zet inundatie = 0 in het insteekvlak
    ds_insteek = rasterio.open(insteek_file)
    ds_inund   = rasterio.open(input_file)
    ds_prep = np.where(ds_insteek.read()[0]==1, 0, ds_inund.read()[0])    
    with rasterio.open(output_file, 'w', **ds_inund.profile) as dest:
        dest.write_band(1,ds_prep)
    ds_inund.close()
    ds_insteek.close()
    
def inundation_total(input_file, output_file, water_level):
    volume_step = add_water_to_raster(input_file, output_file, water_level)
    return volume_step
    
def inundation_excl_allowable(input_file, output_file, topo_file, shp_insteek, insteek_file, extent, water_level):
    # Get shp_insteek
    shp_insteek = gpd.read_file(shp_insteek)
    extent      = gpd.GeoDataFrame([extent])    
    extent.crs  = shp_insteek.crs 
    shp_insteek = gpd.overlay(extent, shp_insteek, how='intersection',keep_geom_type=True)
    shp_insteek['Value'] = 1
    
    # Open height data
    raster     = rasterio.open(input_file)
    meta       = raster.meta.copy()
    meta.update(compress='lzw')
    height_data = raster.read(1)
        
    # Rasterize shp_insteek
    if len(shp_insteek)>0:
        
        with rasterio.open(insteek_file, 'w+', **meta) as out:
            out_arr = out.read(1)
            shapes = ((geom,value) for geom, value in zip(shp_insteek.geometry, shp_insteek.Value))
            burned = rasterio.features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
            out.write_band(1, burned) 
        
        # Read the raster data as a numpy array
        raster_instk= rasterio.open(insteek_file)
        
        # Edit height data (remove area below 10%-percentile)
        cond1           = height_data>-900
        cond2           = raster_instk.read(1)!=1
        height_filtered = height_data[np.where(cond1 & cond2)]
        height_10p      = np.percentile(height_filtered,10)   
        height_filled   = np.where(height_data < height_10p, height_10p, height_data)
        height_filled   = np.where(height_data < -900, height_data, height_filled) # no-data pixels need to remain no-data pixels
        
        with rasterio.open(topo_file, 'w', **raster.profile) as dest:
            dest.write_band(1,height_filled)
        
        raster_instk.close()
    
    else:
        with rasterio.open(topo_file, 'w', **raster.profile) as dest:
            dest.write_band(1,height_data)
            
    # Estimate inundation
    volume_step = add_water_to_raster(topo_file, output_file, water_level)
    
    # Close rasters
    raster.close()
    
    return volume_step
    
    
def fun(polygons):
    # output paths
    output_folder    = config.GIS_output
    excel_file       = os.path.join(config.output, "MM_INNUN.txt")
    raster_topo      = os.path.join(output_folder, "AHN3_clip_peilgebied.tif")
    raster_inun      = os.path.join(output_folder, "inundation.tif")
    raster_insteek   = os.path.join(output_folder, "insteek.tif")
    shp_insteek_temp = os.path.join(output_folder, "insteek.shp")
    shp_insteek      = config.shp_insteek
    
    # Initialize an empty DataFrame to store the results
    volume_total = pd.DataFrame(columns=['polder_name', 'water_level', 'volume_step', 'mm'])
    
    # Open the raster file
    with rasterio.open(config.raster_prep) as src:
        
        # Read the shapefile using geopandas
        # polygons = gpd.read_file(config.peil_clip)[0:2]
        # polygons = polygons[0:5]
        
        # Calculate polygon areas in square meters
        polygons['area_m2'] = polygons['geometry'].area
        
        # Initialize an empty DataFrame to store the results
        volume_total = pd.DataFrame(columns=['code', 'water_level', 'volume_tot_m3', 'volume_inun_m3','volume_tot_mm','volume_inun_mm'])
        
        # Iterate over each polygon
        for idx, polygon in tqdm(polygons.iterrows()):
            # print(polygons['CODE'][idx])
                        
            # Extract the geometry of the polygon
            geometry    = polygon.geometry
            
            # Extract water level
            water_level = polygon['WLmax']
            
            # Clip the raster with the polygon geometry
            topo_clip, out_transform = mask(src, [geometry], crop=True)
            
            # Update the metadata
            topo_meta = src.meta.copy()
            topo_meta.update({
                "driver": "GTiff",
                "height": topo_clip.shape[1],
                "width": topo_clip.shape[2],
                "transform": out_transform,
            })
            
            # Write the clipped raster to a new file
            with rasterio.open(raster_topo.replace('.tif',f"_{idx}.tif"), "w", **topo_meta) as dest:
                dest.write(topo_clip)
            
            # Calculate total volume
            input_file  = raster_topo.replace('.tif',f"_{idx}.tif")
            output_file = raster_inun.replace('.tif',f"_tot{idx}.tif")
            volume_tot  = inundation_total(input_file, output_file, water_level)            
            
            # Calculate volume (after removing allowable inundation)
            input_file  = raster_topo.replace('.tif',f"_{idx}.tif")
            output_file = raster_inun.replace('.tif',f"_inun{idx}.tif")
            insteek_file= raster_insteek.replace('.tif',f"_{idx}.tif")
            topo_file   = raster_topo.replace('.tif',f"_filled{idx}.tif")
            volume_inun = inundation_excl_allowable(input_file, output_file, topo_file, shp_insteek, insteek_file, polygon, water_level)
            
            if volume_tot is not None:
                # Calculate mm only if volume_step is not None
                volume_tot_mm = round(volume_tot / polygon['area_m2'] * 1000, 1)
                volume_inun_mm= round(volume_inun / polygon['area_m2'] * 1000, 1)
                
                # Append to the DataFrame
                volume_total = pd.concat([volume_total, pd.DataFrame([{
                    'code': polygons['CODE'][idx],
                    'water_level': water_level,
                    'volume_tot_m3': np.round(volume_tot,2),
                    'volume_inun_m3':np.round(volume_inun,2),
                    'volume_tot_mm': np.round(volume_tot_mm,2),
                    'volume_inun_mm': np.round(volume_inun_mm,2),
                }])], ignore_index=True)
            else:
                print(f"Skipping water level {water_level} for {polygons['CODE'][idx]} due to an error in processing.")
            
            # Remove temp files
            file = raster_insteek.replace('.tif',f"_{idx}.tif")
            if os.path.exists(file): os.remove(file)
            file = raster_topo.replace('.tif',f"_{idx}.tif")
            if os.path.exists(file): os.remove(file)
            file = raster_topo.replace('.tif',f"_filled{idx}.tif")
            if os.path.exists(file): os.remove(file)
            
    # Use the to_excel method to export the DataFrame to Excel
    volume_total.to_csv(excel_file, index=False)  # Set index to False if you don't want to export the row numbers as a column
    