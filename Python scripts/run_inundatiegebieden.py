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
from osgeo import gdal, ogr, osr
from rasterio import features
from shapely.geometry import shape
from osgeo_utils.gdal_polygonize import gdal_polygonize
from tqdm import tqdm

class paths():
    # Input
    root = os.getcwd()
    tif_AHN3    = os.path.join(root,'GIS','AHN4_ruw.tif')
    shp_insteek = os.path.join(root,'GIS','BR_VS_Insteekvlak.shp')
    shp_hydro   = os.path.join(root,'GIS','BR_VS_Oppervlaktewater.shp')
    shp_peil    = os.path.join(root,'GIS','BR_VS_Peilgebieden.shp')
    
    # temporary
    tif_inundatie = os.path.join(root,'GIS','Peil','inundatie.tif')
    shp_droognat  = os.path.join(root, 'GIS','Peil','shp_droognat.shp')    
    shp_droognat2 = os.path.join(root, 'GIS','Peil','shp_droognat_clipped.shp') 

    # Output
    shp_nat       = os.path.join(root, 'GIS','Peil','shp_nat.shp')    
    

    
if __name__ == "__main__":
    print('Start')
    
    # Input data
    AHN         = rasterio.open(paths.tif_AHN3)    
    shp_insteek = gpd.read_file(paths.shp_insteek)
    # shp_hydro   = gpd.read_file(paths.shp_hydro)
    shp_extent  = gpd.read_file(paths.shp_peil)
    
    # Peil
    peil = [-2.05, -2.00, -1.95]    
    
    for p_iter in range(0,len(peil)):
        # inundatie
        print(p_iter, 'Stap 1: inundatie')
        inundatie = np.where(AHN.read()<peil[p_iter], 1, 0)
        inundatie = np.where(AHN.read()<-1e10, -999, inundatie)
        profile   = AHN.profile
        with rasterio.open(paths.tif_inundatie.replace('Peil','Peil'+str(p_iter)), 'w', **profile) as dest:
                dest.write(inundatie)
        
        # polygonize
        print(p_iter, 'Stap 2: polygonize')
        raster = gdal.Open(paths.tif_inundatie.replace('Peil','Peil'+str(p_iter)))
        band = raster.GetRasterBand(1)
        drv = ogr.GetDriverByName('ESRI Shapefile')
        outfile = drv.CreateDataSource(paths.shp_droognat.replace('Peil','Peil'+str(p_iter)))
        outlayer = outfile.CreateLayer('polygonized raster', srs = None )
        newField = ogr.FieldDefn('Inundatie', ogr.OFTReal)
        outlayer.CreateField(newField)
        gdal.Polygonize(band, None, outlayer, 0, [])
        outfile = None
        
        # clip
        print(p_iter, 'Stap 3: clip & bereken area')
        shp      = gpd.read_file(paths.shp_droognat.replace('Peil','Peil'+str(p_iter)))
        shp_clip = gpd.clip(shp, shp_extent)
        shp_clip.crs='EPSG:28992'
        shp_clip['area'] = shp_clip.area
        shp_clip.to_file(paths.shp_droognat2.replace('Peil','Peil'+str(p_iter)))
            
        # assign value based on shp_insteek
        print(p_iter, 'Stap 4: verwijder gebieden niet verbonden aan de watergangen')
        shp       = gpd.read_file(paths.shp_droognat2.replace('Peil','Peil'+str(p_iter)))
        intersect = np.zeros(len(shp))    
        for index, orig in tqdm(shp.iterrows(),str(len(shp))):               
            if (orig['Inundatie']==-999) & (orig['area'] > 2e5): continue # handmatige grens...
            if orig['Inundatie']==0: continue # droge gebieden skippen
            if (orig['area'] > 1e4): print(index, orig['area'])    
            count = orig['geometry'].intersects(shp_insteek['geometry']).sum()
            if count > 0:
                intersect[index] = 1
        shp['intersect'] = intersect
        shp = shp.loc[shp['intersect']==1]
        shp = shp.loc[shp['Inundatie']!=0]
        shp.to_file(paths.shp_nat.replace('Peil','Peil'+str(p_iter)))
    
  
  
   