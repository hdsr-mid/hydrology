# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 09:37:42 2024

@author: PetraH

Met deze script worden Sobek resultaten (die al gekoppeld zijn aan een shapefile) gevisualiseerd.
- Waterniveau: punten/nodes
- Debiet: lijnen/lines

Output: 2 png figuren


"""

import os
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors


# new colormap
colors = [(1, 1, 1, 1), (0, 0, 0.5, 1)] # first color is white, last is blue
new_cmap = LinearSegmentedColormap.from_list("Custom", colors, N=20)
colors1  = new_cmap(np.linspace(0, 1, 30))
colors2  = plt.cm.jet(np.linspace(0, 1, 128))
colors   = np.vstack((colors1, colors2))
mymap    = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

class paths():
    root            = os.getcwd()
    shp_node        = os.path.join(root,"Data",'Sobek_WL.shp')
    shp_reach       = os.path.join(root,"Data","Sobek_Q.shp")
    shp_peil        = os.path.join(root,'Data','BR_VS_Peilgebieden.shp')
    
def fun_plot_node():
    # load data
    shp_node      = gpd.read_file(paths.shp_node)
    shp_peil      = gpd.read_file(paths.shp_peil)
    
    # Create figure
    plt.figure(figsize=(10,5))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    ax = plt.subplot2grid((1, 1), (0, 0))
    shp_peil.plot(ax=ax, color='white', edgecolor='black')    
    shp_node.plot(ax=ax, column = 'WL', vmin = 0, vmax = 1, cmap = mymap, legend=True, legend_kwds={"label": "Waterniveau verschil [m]"})
    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    plt.savefig('figs_nodes.png')
    plt.close()
    
    
def fun_plot_line():
    # Input data
    shp_peil      = gpd.read_file(paths.shp_peil)
    shp_reach     = gpd.read_file(paths.shp_reach)
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


if __name__ == "__main__":
    fun_plot_node()
    fun_plot_line()
    