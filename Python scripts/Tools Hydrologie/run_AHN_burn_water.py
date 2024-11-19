# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 09:02:03 2023

@author: PetraH

In deze script worden de volgende stappen gevolgd voor AHN data-preprocessing:
    - clip op gebied (shp_extent)
    - watervlaktes in AHN identificeren en de waarden voor deze pixels gelijk zetten aan het streefpeil

Output:
    Geclipt AHN raster 
    Raster voor watervlaktes met streefpeilen
    Geclipt AHN raster incl streefpeil waarden voor watervlaktes


"""

import os
import numpy as np
import rasterio
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import mapping
import warnings
warnings.filterwarnings("ignore")

class paths():
    ''' define all paths here '''
    # Current folder location
    root   = os.getcwd()

    # Paths: Input files
    raster_file = os.path.join(r'H:\DATA\Service\Petra\00 GIS DATA', "AHN4 ruw compressed.tif")
    shp_extent  = os.path.join(root,'Input','EvS_gebied.shp')
    shp_water   = os.path.join(root,'Input','BR_VS_Watervlak.shp')
    shp_peil    = os.path.join(root,'Input','BR_VS_Peilgebieden.shp')
    
    # Paths: output files
    raster_clip = os.path.join(root,'Output',"AHN4_ruw_clip.tif")
    raster_peil = os.path.join(root,'Output','Peil.tif')
    raster_prep = os.path.join(root,'Output',"AHN4_ruw_clip_prep.tif")
   
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

def fun_clip_raster():
    topo       = rasterio.open(paths.raster_file)
    extent     = gpd.read_file(paths.shp_extent)    
    geometries = extent.geometry.apply(mapping)

    topo_clip, out_transform = mask(topo, geometries, crop=True)
    out_meta = get_updated_meta(topo, topo_clip, out_transform)            
    with rasterio.open(paths.raster_clip, 'w', **out_meta) as dest:
        dest.write(topo_clip)

def fun_preprocess_raster():    
   
    # Get water level height from shp_peil and transfer it to shp_water
    print('Raster preparation: get WL')
    extent      = gpd.read_file(paths.shp_extent)    
    shp_water   = gpd.read_file(paths.shp_water)    
    shp_peil    = gpd.read_file(paths.shp_peil)    
    shp_water   = gpd.overlay(extent, shp_water, how='intersection',keep_geom_type=True)
    shp_peil    = gpd.overlay(extent, shp_peil, how='intersection',keep_geom_type=True)
    shp_peil['PEIL'] = np.where(shp_peil["ZOMERPEIL"]==0,shp_peil["VASTPEIL"],shp_peil["ZOMERPEIL"])
    shp_peil['PEIL'] = np.where(shp_peil["PEIL"]==0,shp_peil["FLEXIBEL_B"],shp_peil["PEIL"])
    shp_peil['PEIL'] = np.where(shp_peil['PEIL'] ==-999,np.nan,shp_peil['PEIL'] )
    shp_water  = shp_water.sjoin(shp_peil, how="left")[['geometry','PEIL']]
     
    # Rasterize water level height    
    print('Raster preparation: rasterize WL')
    ds_topo    = rasterio.open(paths.raster_clip)
    meta       = ds_topo.meta.copy()
    meta.update({"compress":'lzw',"nodata": -999})
    with rasterio.open(paths.raster_peil, 'w+', **meta) as out:
        out_arr = out.read(1)
        shapes = ((geom,value) for geom, value in zip(shp_water.geometry, shp_water.PEIL))
        burned = rasterio.features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)
    ds_topo.close()
    
    # Set water level heights in AHN
    print('Raster preparation: burn WL')
    ds_Peil = rasterio.open(paths.raster_peil)
    ds_topo = rasterio.open(paths.raster_clip)
    ds_prep = np.where(ds_Peil.read()[0]==-999, ds_topo.read()[0],ds_Peil.read()[0])
    meta.update({"nodata": -999})
    with rasterio.open(paths.raster_prep, 'w', **ds_topo.profile) as dest:
        dest.write_band(1,ds_prep)
    ds_topo.close()
    ds_Peil.close()
    
    print('Raster preparation: done')
    
if __name__ == "__main__":
    # Prepare inundationmap: Clip AHN raster to target region    
    print('Clip raster')
    fun_clip_raster()
    
    # Prepare inundationmap: AHN gapfilling & set height water level (i.e., op de watergangen is de hoogte gelijk aan het waterpeil)
    print('Preprocess raster')
    fun_preprocess_raster()