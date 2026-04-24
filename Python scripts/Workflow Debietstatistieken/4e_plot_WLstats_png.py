# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 14/04/2026
virtual environment used: geo-env
"""

import os
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")


class general():
    outdir              = r'D:/workingdir/3_Output'
    shp_watergangen     = r'D:/workingdir/1_InputData/Branches.gpkg'
    path_WLstats        = r'D:/workingdir/3_Output/gpkg/Waterstandstatistieken.gpkg'

def fun_plot_png(gdf,stat):
    watergangen = gpd.read_file(general.shp_watergangen)
    
    if 'dWL' in stat:
        vminmax    = [-1,1]
        thresholds = [0.1, 0.5]
        title      = stat[1:] + ' - WL(t=0)'
    else:
        vminmax    = [-5,5]
        thresholds = [1, 2.5]
        title      = 'WL(T=1)'
    size = np.where(abs(gdf[stat])<thresholds[0],0.1,
                    np.where((abs(gdf[stat])>=thresholds[0]) & (abs(gdf[stat])<thresholds[1]),1,
                    2))
    
    plt.figure(figsize=(6,3))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    ax = plt.subplot2grid((1, 1), (0, 0))
    watergangen.plot(ax=ax, color='grey',linewidth=0.3,alpha=1)
    gdf.plot(ax=ax,column=stat,cmap='jet',vmin=vminmax[0],vmax=vminmax[1], legend=True,markersize=size, legend_kwds={"label": stat + " [m]"},alpha=1)
    plt.title(title + ' = ' + str(np.round(gdf[stat].min(),1)) + ' - ' + str(np.round(gdf[stat].max(),1)) + ' m')
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.grid(linestyle='-', linewidth=0.1)
    plt.savefig(os.path.join(general.outdir,'png','fig_'+stat+'.png'),dpi=300)
    plt.close()


if __name__ == "__main__":
    
    # load data
    print('Loading data...')
    gdf = gpd.read_file(general.path_WLstats)
    
    # plot
    print('Plotting...')
    fun_plot_png(gdf, 'dWLmax')
    fun_plot_png(gdf,'dWLmin')
    fun_plot_png(gdf,'WL_T1')
    
    print('Done!')