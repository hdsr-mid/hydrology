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
    dWL_abs    = abs(gdf_joined['dWL'])
    dWL_norm   = dWL_abs/np.max(dWL_abs)
    msize = np.where(abs(gdf_joined['dWL'])<0.01, 5, 5+dWL_norm*50) 
    
    plt.figure(figsize=(10,10))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.1)
    
    ax = plt.subplot2grid((2, 1), (0, 0))
    shp_afvoer.plot(ax=ax, color='white', edgecolor='black')
    gdf_joined.plot(ax=ax, column = 'dWL', vmin = -1, vmax = 1, markersize=msize, cmap = 'jet', legend=True, legend_kwds={"label": "Verschil: Sobek - FEWS_WIS data [m]"})
    plt.xlim([109000,170000])
    plt.ylim([438000,468000])
    plt.grid()
    
    ax = plt.subplot2grid((2, 1), (1, 0))
    plt.title('Op basis van maximale absolute waarde per afvoergebied...')
    shp_afvoer.plot(ax=ax, column = 'dWL', vmin = -1, vmax = 1, cmap = 'jet', edgecolor='black', linewidth=0.1, alpha = 0.5, legend=True, legend_kwds={"label": "Maximum absoluut verschil: Sobek - HYDROMEDAH [m3/s]"})
    shp_afvoer.plot(ax=ax, color= "none", facecolor = "none", edgecolor='black')
    plt.xlim([109000,170000])
    plt.ylim([438000,468000])
    plt.grid()
    
    plt.savefig(paths.fig, dpi=300)
    plt.show()
    plt.close()

def bokeh(paths, shp_afvoer, shp_points):  
    metric = 'dWL'
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
    shp_points    = shp_points.rename(columns ={'X_left':'X','Y_left':'Y'})
    
    p1 = figure(title='Sobek vs. FEWS_WIS debiet data [m]: ' + metric, height=350, width=820)
    p1.patches(fill_alpha=0.7,fill_color='white',line_color='black', line_width=0.5, source=geo_source)
    map_source = ColumnDataSource(shp_points)    
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_points = p1.scatter('X', 'Y', source=map_source,color=transform(metric, color), size='msize')
    color_bar = ColorBar(color_mapper=color,title=metric)
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@ID_right'),
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
    
def fun_T1(x):
     x = x[np.isnan(x)==0]
     x = np.sort(x)[::-1]
     N = len(x)
     i = np.arange(1,N+1)
     p = i/(N+1)
     T = 1/p
     
     if N > 3*365:
         interp_func = interp1d(T,x)
         x_T1        = interp_func(365)
     else:
         x_T1 = -999
     
     return np.round(float(x_T1),2)

def fun_stat_data(df_sites, df_data):
    
    # selecteer data
    xcoord = []
    ycoord = []
    df_stat= pd.DataFrame()
    for xy, group in df_sites.groupby(['X','Y']):
        
        # selecteer locaties die inbegrepen zijn in de data csv
        ID = list(group['LOC_ID'].values)
        group = group.set_index('LOC_ID')
        for i in ID:
            if i not in np.array(df_data.columns):
                group = group.drop(i)
        group = group.reset_index()
        ID = list(group['LOC_ID'].values)
        
        # extract data
        if len(ID)>0:
            df_sel = df_data[np.append('time',ID)]
            for i_ID in ID: df_sel = df_sel.astype({i_ID: 'float'})            
            
            for i_ID in ID:
                data        = df_sel[[i_ID]].replace(-999,np.nan).dropna()
                if len(data)> t_min_coverage:
                    data_stat   = float(data.apply(fun_T1).values)
                else:
                    data_stat = -999
                df_group    = pd.DataFrame({'ID': i_ID,'X':group['X'].values, 'Y': group['Y'].values, 'Data_T1': data_stat})
                df_stat     = df_stat._append(df_group)
    
    return df_stat
 
def main(paths, months, years):
    ''''FEWS-WIS data'''
    # Input data
    df_sites = pd.read_csv(paths.df_WL_WISsites)    
    df_data  = pd.read_csv(paths.df_WL_WISdata, delimiter=',',low_memory=False)    
    df_data  = df_data.rename(columns={"Unnamed: 0": "time"}).drop(0)
    df_data  = df_data.astype({'time':'datetime64[ns]'})    
    df_data['MONTH'] = np.array([pd.to_datetime(t).month for t in df_data.time.values])
    df_data['YEAR'] = np.array([pd.to_datetime(t).year for t in df_data.time.values])
    df_data    = df_data[df_data['MONTH'].isin(months)]
    df_data    = df_data[df_data['YEAR'].isin(years)]
    shp_insteek= gpd.read_file(paths.shp_insteek)
    shp_afvoer = gpd.read_file(paths.shp_afvoer)
    
    # bereken statistiek
    data_stat  = fun_stat_data(df_sites, df_data)
    
    # omzetten naar shapefile
    gdf_data  = gpd.GeoDataFrame(data_stat, geometry=[Point(xy) for xy in zip(data_stat.X,data_stat.Y)])
    gdf_data.crs='EPSG:28992'
    
    # clip op primaire watergangen    
    shp_insteek = shp_insteek.loc[shp_insteek['CATEGORIEO']==1]
    gdf_data = gpd.clip(gdf_data, shp_insteek)
    
    # opslaan als shapefile
    gdf_data.to_file(paths.gdf_data)
    
    
    '''Model data'''
    # data inlezen
    shp_node  = gpd.read_file(paths.shp_node)
    df_WL         = pd.read_csv(paths.df_model,skiprows=4, delimiter = ',').reset_index()
    df_WL.columns = ['index','ID','WL']    
    
    # Link Sobek data to shapefile: Waterniveau
    nmax         = len(shp_node['ID'].values)
    shp_node_WL  = np.zeros((nmax, 1))-999
    for i in range(0,len(df_WL)):
        ind = np.where((df_WL['ID'][i]==shp_node['ID'].values))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_node_WL[ind]  = df_WL['WL'].values[i]         
    shp_node['Model_T1']  = shp_node_WL
    shp_node = shp_node[shp_node['Model_T1'].notna()]
    shp_node.to_file(paths.gdf_model)    
    
    '''Data vs. Model'''
    gdf_data = gpd.read_file(paths.gdf_data)
    gdf_model = gpd.read_file(paths.gdf_model)
    
    gdf_joined = gpd.sjoin_nearest(gdf_data, gdf_model, distance_col="distances")
    gdf_joined['dWL'] = np.where(gdf_joined.distances < node_min_distance, gdf_joined['Model_T1'] - gdf_joined['Data_T1'], -999)    
    gdf_joined['dWL'] = np.where((gdf_joined['Model_T1'] == -999) | (gdf_joined['Data_T1'] == -999), -999, gdf_joined['dWL'])        
    gdf_joined        = gdf_joined[gdf_joined['dWL']!=-999]
    gdf_joined.to_file(paths.gdf_joined)
    
    # spatial join (nodes & shapefile)
    shp_node_joined = gpd.sjoin(gdf_joined[['dWL','geometry']], shp_afvoer[['CODE','geometry']], how="left")
    shp_node_joined['dWLabs'] = abs(shp_node_joined['dWL'])
    shp_node_joined = shp_node_joined.loc[shp_node_joined.groupby('CODE')['dWLabs'].idxmax().values]
    shp_afvoer      = shp_afvoer.merge(shp_node_joined[['CODE','dWL']], on='CODE', how='left')
    
    # Data plotten
    png(paths,gdf_joined,shp_afvoer)
    bokeh(paths, shp_afvoer, gdf_joined)
    