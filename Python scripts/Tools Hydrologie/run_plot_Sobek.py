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
import warnings
warnings.filterwarnings("ignore")
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.patches import Patch

# new colormap
colors = [(1, 1, 1, 1), (0, 0, 0.5, 1)] # first color is white, last is blue
new_cmap = LinearSegmentedColormap.from_list("Custom", colors, N=20)
colors1  = new_cmap(np.linspace(0, 1, 30))
colors2  = plt.cm.jet(np.linspace(0, 1, 128))
colors   = np.vstack((colors1, colors2))
mymap    = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

class paths():
    root            = os.getcwd()
    shp_resultaten  = os.path.join(root,"Output",'resultaten.gpkg')
    lyr_node        = 'Sobek_WL'
    lyr_reach       = 'Sobek_Q'
    lyr_peil        = 'Peilgebieden'
    shp_insteek     = os.path.join(root,'Input','BR_VS_Insteekvlak.shp')  
    
    
def fun_plot_node():
    # load data
    shp_node      = gpd.read_file(paths.shp_resultaten,layer=paths.lyr_node)
    shp_peil      = gpd.read_file(paths.shp_resultaten,layer=paths.lyr_peil)
    
    # Create figure
    plt.figure(figsize=(10,5))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    ax = plt.subplot2grid((1, 1), (0, 0))
    shp_peil.plot(ax=ax, color='white', edgecolor='black')    
    shp_node.plot(ax=ax, column = 'WL', vmin = 0, vmax = 2, cmap = mymap, legend=True, legend_kwds={"label": "Waterniveau [m NAP]"})
    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    plt.savefig('figs_nodes.png')
    plt.close()
        
def fun_plot_line():
    # Input data
    shp_peil      = gpd.read_file(paths.shp_resultaten, layer=paths.lyr_peil)
    shp_reach     = gpd.read_file(paths.shp_resultaten, layer=paths.lyr_reach)
    Q_upper       = 1
    Q_lower       = 0
    
    # Create figure
    plt.figure(figsize=(10,5))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    ax = plt.subplot2grid((1, 1), (0, 0))
    shp_peil.plot(ax=ax, color='white', edgecolor='black')
    ax = shp_reach.plot(ax=ax, column = 'Q', vmin = Q_lower, vmax = Q_upper, cmap = 'jet',linewidth=shp_reach['lw'], legend = True, legend_kwds={"label": "Q [m$^3$/s]"})
    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    
    plt.savefig('fig_line.png')
    plt.close()

def fun_discrete_colormap(data_plotted):
    # bin plotted data
    data_plotted[data_plotted < 2] = 0.1
    data_plotted[(data_plotted < 4) & (data_plotted >= 2)] = 0.2
    data_plotted[(data_plotted < 6) & (data_plotted >= 4)] = 0.3
    data_plotted[(data_plotted < 8) & (data_plotted >= 6)] = 0.4
    data_plotted[(data_plotted < 10) & (data_plotted >= 8)] = 0.5
    data_plotted[(data_plotted >= 10)] = 0.6
    
    # Let's also design our color mapping: 1s should be plotted in blue, 2s in red, etc...
    col_dict = {0.1: "blue", 0.2: "red", 0.3: "orange", 0.4: "green", 0.5: "green", 0.6: "green"}        
    # We create a colormar from our list of colors
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])        
    # Let's also define the description of each category
    labels = np.array(["< 2 cm", "2 - 4 cm", '4 - 6 cm', "6 - 8 cm", "8 - 10 cm",'> 10 cm'])
    len_lab = len(labels)
    # prepare normalizer
    # Prepare bins for the normalizer
    norm_bins = np.sort([*col_dict.keys()]) + 0.5
    norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
    
    # Make normalizer and formatter
    norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
    
    return data_plotted, norm, fmt

def fun_plot_polygon():
    # input data
    shp_peil = gpd.read_file(paths.shp_resultaten, layer=paths.lyr_peil)                
    
    # Create figure
    fig = plt.figure(figsize=(6,2.5))
    plt.subplots_adjust(left=0.01, bottom=0.1, right=0.8, top=0.9, wspace=0.01, hspace=0.9)
    ax = plt.subplot2grid((1, 1), (0, 0))
    
    thresholds = [-1e10, 0, 0.5, 1, 2, 3, 1e10]        
    colors     = ["lightgrey", "navy", "dodgerblue", "gold", "darkorange", "darkred"]
    labels     = np.array(["<0 m NAP", "0 - 0.5 m NAP", '0.5 - 1.0 m NAP', "1.0 - 2.0 m NAP", "2.0 - 3.0 m NAP",'> 2 m NAP'])
    pmarks     = []
    shp_peil.plot(ax=ax, color='white', edgecolor='black', linewidth=0.5)
    for ii in range(0,len(thresholds)-1):
        data_plotted = np.where((shp_peil['WLmax'] >= thresholds[ii]) & (shp_peil['WLmax'] < thresholds[ii+1]), ii, np.nan)
        shp_peil['data_plotted'] = data_plotted
        im = shp_peil.dropna(subset = ['data_plotted']).plot(ax=ax, column = 'data_plotted', edgecolor='black', linewidth=0.5, color = colors[ii], alpha=0.7)    
        pmarks.append(Patch(facecolor=colors[ii], label=labels[ii], alpha=0.7))

    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    plt.xticks([])
    plt.yticks([])
    
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(title = 'Max waterniveau [m NAP]',handles=[*handles,*pmarks], bbox_to_anchor=(1.5, 1.05), title_fontsize=8, fontsize=8)
    
    plt.savefig('fig_polygon.png', dpi=300)
    plt.close()
        
if __name__ == "__main__":
    fun_plot_node()
    fun_plot_line()
    fun_plot_polygon()
    