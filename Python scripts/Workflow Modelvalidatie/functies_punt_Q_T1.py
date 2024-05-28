# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 13:14:34 2024

@author: PetraH
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from tqdm import tqdm
import matplotlib.dates as mdates
from datetime import datetime
import xarray as xr
import warnings
from shapely.geometry import Point
warnings.filterwarnings("ignore")

from bokeh.plotting import figure, save
from bokeh.models import ColumnDataSource, GeoJSONDataSource, LinearColorMapper, ColorBar
from bokeh.layouts import layout
from bokeh.transform import transform
from bokeh.models import HoverTool 

node_min_distance = 5 # min distance between model node and observation point [m]
t_min_coverage    = 3*365 # minimum number observation points required to validate

def png(paths,gdf_joined,shp_afvoer):
    dQ_abs    = abs(gdf_joined['dQ'])
    dQ_norm   = dQ_abs/np.max(dQ_abs)
    msize     = np.where(abs(gdf_joined['dQ'])<0.01, 5,5 + dQ_norm*50) 
    
    
    plt.figure(figsize=(10,10))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.1)
    
    ax = plt.subplot2grid((2, 1), (0, 0))
    shp_afvoer.plot(ax=ax, color='white', edgecolor='black')
    gdf_joined.plot(ax=ax, column = 'dQ', vmin = -1, vmax = 1, markersize=msize, cmap = 'jet', legend=True, legend_kwds={"label": "Verschil: Sobek - FEWS_WIS data [m3/s]"})
    plt.xlim([109000,170000])
    plt.ylim([438000,468000])
    plt.grid()
    
    ax = plt.subplot2grid((2, 1), (1, 0))
    plt.title('Op basis van maximale absolute waarde per afvoergebied...')
    shp_afvoer.plot(ax=ax, column = 'dQ', vmin = -1, vmax = 1, cmap = 'jet', edgecolor='black', linewidth=0.1, alpha = 0.5, legend=True, legend_kwds={"label": "Maximum absoluut verschil: Sobek - HYDROMEDAH [m3/s]"})
    shp_afvoer.plot(ax=ax, color= "none", facecolor = "none", edgecolor='black')
    plt.xlim([109000,170000])
    plt.ylim([438000,468000])
    plt.grid()
    
    plt.savefig(paths.fig, dpi=300)
    plt.show()
    plt.close()

def bokeh(paths, shp_afvoer, shp_points):  
    metric = 'dQ'
    vmin   = -1
    vmax   = 1
    
    # -----------------------------------------------
    # Get data
    shp_points = shp_points[shp_points[metric].replace(np.nan,-999)!=-999]
    shp_points = shp_points.drop('geometry', axis=1).copy()
    geo_source = GeoJSONDataSource(geojson=shp_afvoer.to_json())
    
    # -----------------------------------------------
    # Plot map: with points
    metric_abs    = abs(shp_points[metric])
    metric_norm   = metric_abs/np.max(metric_abs)
    shp_points['msize'] = np.where(abs(shp_points[metric])<0.05, 0.05, 1+metric_norm *20)     
    
    p1 = figure(title='Sobek vs. FEWS_WIS debiet data [m3/s]: ' + metric, height=350, width=820)
    p1.patches(fill_alpha=0.7,fill_color='white',line_color='black', line_width=0.5, source=geo_source)
    map_source = ColumnDataSource(shp_points)    
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_points = p1.scatter('X', 'Y', source=map_source,color=transform(metric, color), size='msize')
    color_bar = ColorBar(color_mapper=color,title=metric)
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@IRIS_ID'),
                ('Label', '@'+metric)]
    hover = HoverTool(renderers=[map_points], tooltips=tooltips) 
    p1.add_tools(hover)  
    
    # -----------------------------------------------
    # Plot map: without points
    p2 = figure(title='Op basis van maximale absolute waarde per afvoergebied...', height=350, width=700)
    map_poly = p2.patches(fill_alpha=0.7,
              fill_color={'field': metric, 'transform': color},
              line_color='black', line_width=0.5, source=geo_source)
    tooltips = [('Code', '@NAAM'),
                ('Label', '@'+metric)]
    hover = HoverTool(renderers=[map_poly], tooltips=tooltips) 
    p2.add_tools(hover)  
    
    grid = layout([[p1, p2],])
    save(grid,paths.fightml)
    
def fun_remove(df, name):
    df_sel = df.copy()
    for xy, group in df.groupby(['X','Y']):
        if (name in group['TYPE'].values) & (len(group['TYPE'].values)>0):
            ind = np.where((df_sel['X']==xy[0]) & (df_sel['Y']==xy[1]) & (df_sel['TYPE']!=name))[0]
            ind = df_sel.iloc[ind].index.values            
            df_sel = df_sel.drop(ind)
            
    return df_sel

def add_stats(df_sites_Q, df_data_Q):
    
    # selecteer data
    for i in range(0,len(df_sites_Q)):            
        # selecteer locaties die inbegrepen zijn in de data csv
        ID = df_sites_Q.iloc[i]['ID']
        
        # extract data & add to xds
        df_sel = df_data_Q[np.append('time',ID)]
        df_sel = df_sel.astype({ID: 'float'})                               
         
        # Estimate stats
        data_T1     = fun_T1(df_sel[ID].replace(-999,np.nan))  
        
        # add to dataframe
        df_sites_Q.loc[(df_sites_Q['ID'] == ID), 'Q_T1_data'] = data_T1
    
    return df_sites_Q

def fun_T1(x):
    x = x[np.isnan(x)==0]
    x = np.sort(x)[::-1]
    N = len(x)
    i = np.arange(1,N+1)
    p = i/(N+1)
    T = 1/p
    
    if N > t_min_coverage:
        interp_func = interp1d(T,x)
        x_T1        = interp_func(365)
    else:
        x_T1 = -999
    
    return np.round(float(x_T1),2)

def Sobek_2_shp_Q(df_Q, shp_reach):
    # Link Sobek data to shapefile: Debiet
    rmax        = len(shp_reach['ID'].values)
    shp_reach_Q = np.zeros((rmax,))*np.nan
    for i in range(0,len(shp_reach_Q)):
        ind = np.where((shp_reach['ID'].values[i]==df_Q['ID']))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_reach_Q[i] = float(df_Q['Q_T1_Sobek'].values[ind])
    shp_reach['Q_T1_Sobek'] = shp_reach_Q
    shp_reach = shp_reach[shp_reach['Q_T1_Sobek'].notna()]
    
    return shp_reach

def fun_remove_duplicates(gdf_data):
    # Step 1 - collect all rows that are *not* duplicates (based on ID)
    non_duplicates_to_keep = gdf_data.drop_duplicates(subset='ID_left', keep=False)

    # Step 2a - identify *all* rows that have duplicates (based on ID, keep all)
    sub_df = gdf_data[gdf_data.duplicated('ID_left', keep=False)]

    # Step 2b - of those duplicates, discard all that have "0" in column
    duplicates_to_keep = sub_df[sub_df['Q_T1_Sobek'] != 0]

    # join the 2 sets
    gdf_data = pd.concat([non_duplicates_to_keep, duplicates_to_keep])
    
    # sort
    gdf_data = gdf_data.sort_values(by=['ID_left'], ascending=True)
    
    return gdf_data
    
def main(paths, months, years):
    ''''FEWS-WIS data'''
    # Input data
    df_sites_Q = pd.read_csv(paths.df_Q_WISsites)   
    df_sites_Q = df_sites_Q.rename(columns={"LOC_ID": "ID"}).drop(0)
    df_data_Q  = pd.read_csv(paths.df_Q_WISdata, delimiter=',',low_memory=False)
    df_data_Q  = df_data_Q.rename(columns={"Unnamed: 0": "time"}).drop(0)
    df_data_Q  = df_data_Q.astype({'time':'datetime64[ns]'}) 
    df_data_Q['MONTH'] = np.array([pd.to_datetime(t).month for t in df_data_Q.time.values])
    df_data_Q['YEAR'] = np.array([pd.to_datetime(t).year for t in df_data_Q.time.values])
    df_data_Q  = df_data_Q[df_data_Q['MONTH'].isin(months)]
    df_data_Q  = df_data_Q[df_data_Q['YEAR'].isin(years)]
    shp_insteek= gpd.read_file(paths.shp_insteek)
    shp_afvoer = gpd.read_file(paths.shp_afvoer)
    
    # Select type: debietmeter, pompvijzel
    df_sites_Q = df_sites_Q.loc[(df_sites_Q['TYPE']=='debietmeter') | (df_sites_Q['TYPE']=='pompvijzel')]
    idx_data   = [i for i in range(0,len(df_sites_Q)) if df_sites_Q.ID.values[i] in df_data_Q.columns]
    df_sites_Q = df_sites_Q.iloc[idx_data] 
    
    # Get data
    data_T1  = add_stats(df_sites_Q, df_data_Q) # selecteer data
    
    # Create geodataframe
    gdf_data  = gpd.GeoDataFrame(data_T1, geometry=[Point(xy) for xy in zip(data_T1.X,data_T1.Y)])
    gdf_data.crs='EPSG:28992'
    
    # Cclip op primary waterways
    shp_insteek = shp_insteek.loc[shp_insteek['CATEGORIEO']==1]
    gdf_data    = gpd.clip(gdf_data, shp_insteek)
    
    # Save data
    gdf_data.to_file(paths.gdf_data)
    
    '''Model data'''
    # data inlezen
    shp_reach    = gpd.read_file(paths.shp_reach)
    df_Q_model   = pd.read_csv(paths.df_Q_model,skiprows=4, sep = ',')
    df_Q_model.columns = ['ID','Q_T1_Sobek']
    df_Q_model   = df_Q_model.astype({'ID':'str','Q_T1_Sobek':'float'}) 
    
    # Spaties in kolom namen verwijderen
    shp_reach.columns = [col.strip() for col in shp_reach.columns]
    
    # Link Sobek data aan shapefile
    gdf_model = Sobek_2_shp_Q(df_Q_model, shp_reach)
    
    # Save file
    gdf_model.crs='EPSG:28992'
    gdf_model.to_file(paths.gdf_model)
    
    '''Data vs. Model'''
    gdf_data = gpd.read_file(paths.gdf_data)
    gdf_model = gpd.read_file(paths.gdf_model)
    
    gdf_joined = gpd.sjoin_nearest(gdf_data, gdf_model, distance_col="distances")
    gdf_joined['dQ'] = gdf_joined['Q_T1_Sobek'] - gdf_joined['Q_T1_data']
    # gdf_joined['dQ'] = np.where(gdf_joined.distances < node_min_distance, gdf_joined['dQ'], -999)    
    gdf_joined['dQ'] = np.where((gdf_joined['Q_T1_Sobek'] == -999) | (gdf_joined['Q_T1_data'] == -999), -999, gdf_joined['dQ'])        
    gdf_joined        = gdf_joined[gdf_joined['dQ']!=-999]
    
    # Remove duplicates (remove row with Q=0)
    gdf_joined = fun_remove_duplicates(gdf_joined)
    
    # Save file
    gdf_joined.to_file(paths.gdf_joined)
    
    # spatial join (nodes & shapefile)
    shp_node_joined = gpd.sjoin(gdf_joined[['dQ','geometry']], shp_afvoer[['CODE','geometry']], how="left")
    shp_node_joined['dQabs'] = abs(shp_node_joined['dQ'])
    shp_node_joined = shp_node_joined.loc[shp_node_joined.groupby('CODE')['dQabs'].idxmax().values]
    shp_afvoer      = shp_afvoer.merge(shp_node_joined[['CODE','dQ']], on='CODE', how='left')
    
    # Data plotten
    png(paths,gdf_joined,shp_afvoer)
    bokeh(paths, shp_afvoer, gdf_joined)
    

    
    
    
    