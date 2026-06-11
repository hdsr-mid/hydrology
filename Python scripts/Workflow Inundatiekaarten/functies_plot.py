# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 09:37:42 2024

@author: PetraH
"""

import os
import numpy as np
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")
import config
import rasterio
from rasterio.plot import show
from rasterstats import zonal_stats
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.patches import Patch

CODE = config.CODE

# new colormap
colors = [(1, 1, 1, 1), (0, 0, 0.5, 1)] # first color is white, last is blue
new_cmap = LinearSegmentedColormap.from_list("Custom", colors, N=20)
colors1  = new_cmap(np.linspace(0, 1, 30))
colors2  = plt.cm.jet(np.linspace(0, 1, 128))
colors   = np.vstack((colors1, colors2))
mymap    = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
          
def fun_plot_inundation(shp_peil, raster_file, fig_file):  
    
    # load data
    inundation       = rasterio.open(raster_file)
    shp_insteek      = gpd.read_file(config.shp_insteek)
    # data = np.round(inundation.read()[0]*100,1) # [m] -> [cm]
    data =  inundation.read()[0].astype('float32')
    data[data<=0] = np.nan
    
    # Create figure
    fig=plt.figure(figsize=(10,7))
    # plt.subplots_adjust(left=0.1, bottom=0.01, right=0.9, top=0.9, wspace=0.1, hspace=0.2)
    ax = plt.subplot2grid((1, 1), (0, 0))
    cmax = 40
    image_hidden = ax.imshow(data, cmap=mymap, vmin=0, vmax=cmax)
    show(data, transform=inundation.transform, cmap=mymap, vmin = 0, vmax = cmax, ax=ax, interpolation = 'nearest')
    shp_peil.boundary.plot(ax=ax, edgecolor='black',linewidth=0.5)
    shp_insteek.plot(ax=ax, color='white', linewidth=0)
    cbar = plt.colorbar(image_hidden, ax=ax)
    cbar.set_label('Inundatie [cm]')
    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.savefig(os.path.join(config.output,fig_file), dpi=300)
    plt.close()
    inundation.close()
    
    resolution  = abs(inundation.transform[0] * inundation.transform[4])
    A_tot       = gpd.read_file(config.shp_extent)['geometry'].area.values
    A_inun      = len(np.where(data>0)[0]) * resolution
    data2       = np.where(data > 0, data, 0)  # Only consider positive inundation
    Volume      = np.sum(data2/100) * resolution    
    stat_perc   = np.nanpercentile(data,np.arange(100,50,-5))
    stat_mean   = Volume/A_tot*100
    
    # print('Percentage inundatie oppervlak: ', 100*A_inun/A_tot)
    # print('Inundatie volume: ', Volume)
    # print('Inundatiediepte percentielen: ', stat_perc)
    # print('Inundatiediepte gebiedsgemiddeld: ', stat_mean)
    
    
    stats     = ['percentile_95','mean','max']        
    gdf       = gpd.read_file(config.shp_extent)[[CODE, 'geometry']]
    output_file = 'temp.tif'
    with rasterio.open(raster_file) as dataset:
        # Read the raster data as a numpy array
        data = dataset.read(1)

        # Create a new raster file with inundation values
        with rasterio.open(
            output_file,
            'w',
            driver='GTiff',
            height=dataset.height,
            width=dataset.width,
            count=1,
            compress='lzw',
            dtype='int16',
            crs=dataset.crs,
            transform=dataset.transform,            
        ) as output_dataset:
            output_dataset.write(data, 1)
            
    stat      = zonal_stats(gdf,output_file,stats=stats)
    # print('Inundatiediepte: ', stat)


def fun_plot_inundatie_gem(shp_peil):
    # input data
    inund    = pd.read_csv(os.path.join(config.output, "MM_INNUN.txt"))
    shp_peil['volume_inun_cm'] = inund['volume_inun_mm']/10
    
    selectie = shp_peil[[CODE,'volume_inun_cm']]
    selectie = selectie[selectie[CODE]=='PG0154']
        
    # Create figure
    fig = plt.figure(figsize=(10,2.5))
    plt.subplots_adjust(left=0.01, bottom=0.1, right=0.8, top=0.9, wspace=0.01, hspace=0.9)
    ax = plt.subplot2grid((1, 1), (0, 0))
        
    thresholds = [-1e10, 0.5, 2.5, 5, 7.5, 10, 1e10]        
    colors     = ["lightgrey", "navy", "dodgerblue", "gold", "darkorange", "darkred"]
    labels     = np.array(["<0.5 cm", "0.5 - 2.5 cm", '2.5 - 5 cm', "5 - 7.5 cm", "7.5 - 10 cm",'> 10 cm'])
    pmarks     = []
    for ii in range(0,len(thresholds)-1):
        data_plotted = np.where((shp_peil['volume_inun_cm'] >= thresholds[ii]) & (shp_peil['volume_inun_cm'] < thresholds[ii+1]), ii, np.nan)
        shp_peil['data_plotted'] = data_plotted
        im = shp_peil.dropna(subset = ['data_plotted']).plot(ax=ax, column = 'data_plotted', edgecolor='black', linewidth=0.5, color = colors[ii], alpha=0.7)    
        pmarks.append(Patch(facecolor=colors[ii], label=labels[ii], alpha=0.7))

    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    plt.xticks([])
    plt.yticks([])
    
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(title = 'Gem. inundatiediepte [cm]',handles=[*handles,*pmarks], bbox_to_anchor=(1.7, 1.05), title_fontsize=8, fontsize=8)
    
    plt.savefig(os.path.join(config.output,'fig_inundatiediepte.png'), dpi=300)
    plt.close()
    pmarks = []
    
