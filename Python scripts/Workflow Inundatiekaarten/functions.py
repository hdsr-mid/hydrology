# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 16:36:40 2023

@author: PetraH
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib.dates as mdates
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")
import config
import functies_inundatie
import functies_WL10p
import functies_plot
import functies_inundatie_p2
import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterstats import zonal_stats
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors

CODE = config.CODE

def defineer_peil(shp_node):
    # Defineer streefpeil
    Peil = np.zeros((len(shp_node),))
    for p in range(0,len(Peil)):
        label = str(shp_node['WS_LABEL'].values[p])
        if ('zp' in label):
            Peil[p] = shp_node['ZOMERPEIL'].values[p]
        elif ('vastpeil' in label):
            Peil[p] = shp_node['VASTPEIL'].values[p]
        elif ('flexibel' in label):
            Peil[p] = shp_node['FLEXIBEL_B'].values[p]
        elif ('geen actief peilbeheer' in label):
            Peil[p] = -999
        else:
            Peil[p] = -999           
    shp_node['Peil'] = Peil
    
    return shp_node

def Sobek_2_shp_WL(df_WL, shp_node):
    # Link Sobek data to shapefile: Waterlevel
    df_WL          = df_WL.drop(['Unnamed: 0','Location:'], axis=1)
    WLmax          = df_WL.max().reset_index()
    WLmax.columns  = ['ID','WL']
    
    nmax           = len(shp_node['ID'].values)
    shp_node_WLmax = np.zeros((nmax, 1))*np.nan    
    for i in range(0,len(shp_node)):
        ind = np.where((shp_node['ID'].values[i]==WLmax['ID']))[0]
        if len(ind)>0: 
            ind               = int(ind)
            shp_node_WLmax[i] = WLmax['WL'].values[ind]            
    shp_node['WLmax']  = shp_node_WLmax    
       
    return shp_node

def point_to_poly(shp_node, shp_poly):
    
    # Add max(WLmax) to shp_poly
    shp_node_max = shp_node.groupby(CODE)['WLmax'].max().reset_index()
    WL_max_peil  = np.zeros(len(shp_poly[CODE]),)-99
    for i in range(0,len(shp_poly[CODE])):
        idx            = np.where(shp_poly[CODE][i]==shp_node_max[CODE])[0]
        if len(idx)>0:
            WL_max_peil[i] = shp_node_max.iloc[idx]['WLmax'].values
    shp_poly['WLmax']  = WL_max_peil
    
    return shp_poly

def main():
    
    # Input data
    shp_peil     = gpd.read_file(config.shp_peil)
    shp_node     = gpd.read_file(config.shp_node)
    df_WL        = pd.read_csv(config.df_WL)
    
    # Edit shapefile column notation
    columns_shp_node  = [c.replace('  ','').replace('  ','').replace(' ','') for c in shp_node.columns]
    shp_node.columns  = columns_shp_node
    
    # Streefpeil
    # shp_node = defineer_peil(shp_node)
    
    # Link Sobek data to shapefile: Waterlevel
    shp_node = Sobek_2_shp_WL(df_WL, shp_node)
    
    # Add threshold (WL at 10% inundation) to shp_node
    shp_poly_10p = functies_WL10p.fun(shp_peil)
    shp_peil['MV_10p']  = shp_poly_10p['MV_10p']
    
    # save file
    shp_node.to_file(config.shp_Sobek_WL)          
    
    # Add WLmax to shp_poly
    shp_poly = point_to_poly(shp_node, shp_peil)
    shp_poly.to_file(config.shp_peil_WLmax)                
    
    # Inundatiekaart maken op basis van maximaal waterniveau
    raster_output1 = config.raster_inun
    raster_output2 = config.raster_inun.replace('.tif','_tot.tif') 
    
    # Create inundationmap
    functies_inundatie.fun(shp_poly)

    # Merge inundation maps
    print('Merge inundation maps')
    tif_files = [os.path.join(config.GIS_output,f) for f in os.listdir(config.GIS_output) if f.startswith("inundation_inun") and f.endswith('tif')]
    inundation_map, out_trans = merge(tif_files)
    src      = rasterio.open(tif_files[0])
    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff","height": inundation_map.shape[1],"width": inundation_map.shape[2],"transform": out_trans,})
    with rasterio.open(raster_output1, "w", **out_meta) as dest:
        dest.write(inundation_map)
    src.close()
    for f in tif_files: os.remove(f)
    
    tif_files = [os.path.join(config.GIS_output,f) for f in os.listdir(config.GIS_output) if f.startswith("inundation") and ('_inun' not in f) and f.endswith('tif')]
    inundation_map, out_trans = merge(tif_files)
    src      = rasterio.open(tif_files[0])
    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff","height": inundation_map.shape[1],"width": inundation_map.shape[2],"transform": out_trans,})
    with rasterio.open(raster_output2, "w", **out_meta) as dest:
        dest.write(inundation_map)
    src.close()
    for f in tif_files: os.remove(f)
      
    # Clean inundation maps
    print('Remove wet areas not connected to waterway')
    functies_inundatie_p2.main(raster_output1)
    functies_inundatie_p2.main(raster_output2)
        
    fig_file1 = 'fig_inundatie.png'
    fig_file2 = 'fig_inundatie_tot.png'
    functies_plot.fun_plot_inundation(shp_peil, raster_output1.replace('.tif','_verbonden.tif'), fig_file1)
    functies_plot.fun_plot_inundation(shp_peil, raster_output2.replace('.tif','_verbonden.tif'), fig_file2)
    functies_plot.fun_plot_inundatie_gem(shp_poly_10p)