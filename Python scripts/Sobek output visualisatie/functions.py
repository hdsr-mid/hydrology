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
import functies_AHN
import functies_inundatie
import rasterio
from rasterio.merge import merge
from rasterio.plot import show

def Sobek_2_shp_Q(t, df_Q, shp_reach):
    # Link Sobek data to shapefile: Debiet
    df_Q        = df_Q.drop(['Unnamed: 0','Location:'], axis=1)
    Q_t         = df_Q.iloc[t].reset_index()
    Q_t.columns = ['ID','Q']
    rmax        = len(shp_reach['ID'].values)
    shp_reach_Q = np.zeros((rmax,))*np.nan
    lwidth      = np.zeros((rmax,))*np.nan
    Q_upper     = np.percentile(df_Q.values,95)
    Q_lower     = np.percentile(df_Q.values,5)
    upper_limit = 5
    for i in range(0,len(shp_reach_Q)):
        ind = np.where((shp_reach['ID'].values[i]==Q_t['ID']))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_reach_Q[i] = float(Q_t['Q'].values[ind])
            l_width        = (float(shp_reach_Q[i]) - Q_lower)/(Q_upper - Q_lower)
            lwidth[i]      = np.nanmin([np.nanmax([0, l_width  * upper_limit]),upper_limit]) + 0.5
    shp_reach['Q'] = shp_reach_Q
    shp_reach['lw']= lwidth
    
    return shp_reach, Q_upper, Q_lower

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
            Peil[p] = shp_node['FLEXIBEL_B'].values[p] ##### check this!!        
        else:
            Peil[p] = np.nan            
    shp_node['Peil'] = Peil
    
    shp_node['awl'] = (shp_node['dl_msl'] - shp_node['Peil'])/2 + shp_node['Peil']
    return shp_node

def Sobek_2_shp_WL(t, df_WL, shp_node, suffix):
    # Link Sobek data to shapefile: Waterlevel
    df_WL          = df_WL.drop(['Unnamed: 0','Location:'], axis=1)
    WL_0           = df_WL.iloc[0].reset_index()
    WL_t           = df_WL.iloc[t].reset_index()
    WLmax          = df_WL.max().reset_index()
    WL_t.columns   = ['ID','WL']
    WL_0.columns   = ['ID','WL']
    WLmax.columns  = ['ID','WL']
    
    nmax           = len(shp_node['ID'].values)
    shp_node_WL    = np.zeros((nmax, 1))*np.nan
    shp_node_WL0   = np.zeros((nmax, 1))*np.nan
    shp_node_WLmax = np.zeros((nmax, 1))*np.nan    
    for i in range(0,len(shp_node_WL)):
        ind = np.where((shp_node['ID'].values[i]==WL_t['ID']))[0]
        if len(ind)>0: 
            ind               = int(ind)
            shp_node_WL[i]    = WL_t['WL'].values[ind]    
            shp_node_WL0[i]   = WL_0['WL'].values[ind]
            shp_node_WLmax[i] = WLmax['WL'].values[ind]            
    shp_node['WL']     = shp_node_WL
    shp_node['WLmax']  = shp_node_WLmax
    shp_node['Rvert']  = shp_node['awl'] -  shp_node['WL']
    if 'init' in suffix:
        shp_node['dWL'] = shp_node['WL'] - shp_node['Peil']            
    else:    
        shp_node['dWL']   = shp_node_WL - shp_node_WL0
        # shp_node['dWL']   = np.where(abs(shp_node['dWL'])<0.05, np.nan, shp_node['dWL'])                
    
    return shp_node

def fun_plot_Sobek(t, time, df_P, shp_reach, shp_node, shp_peil, suffix, df_Q, Q_upper, Q_lower):
    # Create figure
    plt.figure(figsize=(10,15))
    # plt.subplots_adjust(left=0.1, bottom=0.01, right=0.9, top=0.9, wspace=0.1, hspace=0.2)
    ax = plt.subplot2grid((5, 1), (0, 0))
    plt.title(time)
    tijd = [datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in df_P['Tijd'].values]    
    plt.bar(tijd, df_P['P [mm/h]'].values,width=0.1)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.gcf().autofmt_xdate()
    plt.xticks(rotation = 0)
    plt.plot([time,time],[0,10],'r')
    plt.ylabel('P [mm/h]')
    plt.xlabel('Tijd')
    
    ax = plt.subplot2grid((5, 1), (1, 0), rowspan = 2)
    shp_peil.plot(ax=ax, color='white', edgecolor='black')
    ax = shp_reach.plot(ax=ax, column = 'Q', vmin = Q_lower, vmax = Q_upper, cmap = 'jet',linewidth=shp_reach['lw'], legend = True, legend_kwds={"label": "Q [m$^3$/s]"})
    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    
    
    ax = plt.subplot2grid((5, 1), (3, 0), rowspan = 2)
    shp_peil.plot(ax=ax, color='white', edgecolor='black')    
    if suffix == 'init':
        shp_node['dWL']   = np.where(abs(shp_node['dWL'])<0.1, np.nan, shp_node['dWL'])                
        shp_node.plot(ax=ax, column = 'dWL', color = 'grey', markersize =1)
        shp_node.plot(ax=ax, column = 'dWL', vmin = -0.5, vmax = 0.5, cmap = 'jet_r', legend=True, legend_kwds={"label": "Waterniveau verschil [m]"})
    else:        
        shp_node['Rvert'] = np.where(shp_node['Rvert']>0, np.nan, shp_node['Rvert']) 
        shp_node.plot(ax=ax, column = 'Rvert', color = 'grey', markersize =1)
        shp_node.plot(ax=ax, column = 'Rvert', vmin = -0.25, vmax = 0.25, cmap = 'jet_r', legend=True, legend_kwds={"label": "Verticale Ruimte [m]"})
    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    plt.savefig(os.path.join(config.output,'figs','fig_'+suffix + '_t' + str(t)+'.png'))
    plt.close()

def fun_plot_inundation(t, time, df_P, shp_peil, suffix, raster_file, fig_file):    
    inundation       = rasterio.open(raster_file)
    
    # Create figure
    fig=plt.figure(figsize=(10,7))
    # plt.subplots_adjust(left=0.1, bottom=0.01, right=0.9, top=0.9, wspace=0.1, hspace=0.2)
    ax = plt.subplot2grid((1, 1), (0, 0))
    data = np.round(inundation.read()[0]*100,1) # [m] -> [cm]
    data[data==0] = np.nan
    cmax = 75
    image_hidden = ax.imshow(data, cmap='jet', vmin=0, vmax=cmax)
    show(data, transform=inundation.transform, cmap='jet', vmin = 0, vmax = cmax, ax=ax)
    shp_peil.boundary.plot(ax=ax, edgecolor='black',linewidth=0.5)
    cbar = plt.colorbar(image_hidden, ax=ax)
    cbar.set_label('Inundatie [cm]')
    plt.xlim([135000,152500])
    plt.ylim([440000,450000])
    plt.grid()
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.title(suffix)
    
    plt.savefig(os.path.join(config.output,'figs_inundatie',fig_file), dpi=300)
    plt.close()
    inundation.close()

def point_to_poly(shp_node, shp_peil):
    # Add WLmax to shp_peil
    shp_node_max = shp_node.groupby('CODE')['WLmax'].max().reset_index()
    WL_max_peil  = np.zeros(len(shp_peil['CODE']),)-99
    for i in range(0,len(shp_peil['CODE'])):
        idx            = np.where(shp_peil['CODE'][i]==shp_node_max['CODE'])[0]
        if len(idx)>0:
            WL_max_peil[i] = shp_node_max.iloc[idx]['WLmax'].values
    shp_peil['WLmax']  = WL_max_peil
    
    return shp_peil

def main(t,suffix):
       
    # Input data
    shp_peil  = gpd.read_file(config.peil_clip)
    shp_reach = gpd.read_file(config.shp_reach)
    shp_node  = gpd.read_file(config.shp_node)
    df_Q      = pd.read_csv(config.df_Q.replace('.txt','_'+suffix+'.txt'))
    df_WL     = pd.read_csv(config.df_WL.replace('.txt','_'+suffix+'.txt'))
    df_P      = pd.read_csv(config.df_P.replace('.txt','_'+suffix+'.txt'))
    
    # Edit shapefile column notation
    columns_shp_reach = [c.replace('  ','').replace('  ','').replace(' ','') for c in shp_reach.columns]
    columns_shp_node  = [c.replace('  ','').replace('  ','').replace(' ','') for c in shp_node.columns]
    shp_reach.columns = columns_shp_reach
    shp_node.columns  = columns_shp_node
    
    # Get time
    time        = df_Q['Location:'].values[t]    
    # time        = datetime.strptime(time, '%Y/%m/%d %H:%M:%S')
    time        = datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
    
    # Link Sobek data to shapefile: Debiet
    [shp_reach,Q_upper, Q_lower] = Sobek_2_shp_Q(t, df_Q, shp_reach)
    
    # Streefpeil
    shp_node = defineer_peil(shp_node)
    
    # Link Sobek data to shapefile: Waterlevel
    shp_node = Sobek_2_shp_WL(t, df_WL, shp_node, suffix)
    if suffix == 'init': shp_node.to_file(config.Sobek_dWL.replace('.shp','_'+suffix+'.shp'))  
    
    # Add WLmax to shp_peil
    shp_peil = point_to_poly(shp_node, shp_peil)
    
    # Plot
    fun_plot_Sobek(t, time, df_P, shp_reach, shp_node, shp_peil, suffix, df_Q, Q_upper, Q_lower)
       
    # Inundatiekaart maken op basis van maximaal waterniveau
    raster_output1 = config.raster_inun.replace('.tif','_'+ suffix + '.tif')
    raster_output2 = config.raster_inun.replace('.tif','_tot_'+ suffix + '.tif')
    if not os.path.exists(raster_output1):
        # Create inundationmap
        functies_inundatie.fun(shp_peil)

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
    
    if t==0:
        fig_file1 = 'fig_inundatie_'+suffix +'.png'
        fig_file2 = 'fig_inundatie_tot_'+suffix +'.png'
        fun_plot_inundation(t, time, df_P, shp_peil, suffix.replace('_',' '), raster_output1, fig_file1)
        fun_plot_inundation(t, time, df_P, shp_peil, suffix.replace('_',' '), raster_output2, fig_file2)