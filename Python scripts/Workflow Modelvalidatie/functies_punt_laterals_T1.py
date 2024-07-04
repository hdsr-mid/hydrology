# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:20:57 2024

@author: PetraH

gebaseerd op scripts ontvangen voor het berekenen van debietsstatistieken
"""
import os
import numpy as np
from scipy.interpolate import interp1d
import geopandas as gpd
import hkvsobekpy as hkv
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from bokeh.plotting import figure, save, output_file
from bokeh.models import ColumnDataSource, GeoJSONDataSource, LinearColorMapper, ColorBar
from bokeh.layouts import layout
from bokeh.transform import transform
from bokeh.models import HoverTool 

def png(paths,shp_afwatr,shp_afvoer,shp_comp):
    # Smaller marker size for when the difference is insignificant
    dQ_abs    = abs(shp_comp['dQ'])
    dQ_norm   = dQ_abs/np.max(dQ_abs)
    msize = np.where(abs(shp_comp['dQ'])<0.01, 1, 1+dQ_norm*5)    
    
    plt.figure(figsize=(10,5))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    ax = plt.subplot2grid((1, 1), (0, 0), rowspan = 2)
    shp_afwatr.plot(ax=ax, column = 'dQ', vmin = -0.1, vmax = 0.1, markersize=msize, cmap = 'jet', edgecolor='face', linewidth=0, alpha = 1, legend=True, legend_kwds={"label": "Verschil: Sobek - HYDROMEDAH [m3/s]"})
    # shp_comp.plot(ax=ax, column = 'dQ', vmin = -0.1, vmax = 0.1, markersize=msize, cmap = 'jet', legend=True, legend_kwds={"label": "Verschil: Sobek - HYDROMEDAH [m3/s]"})
    shp_afvoer.plot(ax=ax, color= "none", facecolor = "none", edgecolor='black')
    plt.xlim([109000,170000])
    plt.ylim([438000,468000])
    plt.grid()
    plt.savefig(paths.fig, dpi=300)
    plt.show()
    plt.close()

def bokeh(paths, shp_afvoer, shp_afwatr):  
    metric = 'dQ'
    vmin   = -0.1
    vmax   = 0.1
    
    # -----------------------------------------------
    # Get data
    # shp_afwatr   = shp_afwatr[shp_afwatr[metric].replace(np.nan,-999)!=-999]
    # shp_afwatr   = shp_afwatr.drop('geometry', axis=1).copy()
    geo_source = GeoJSONDataSource(geojson=shp_afvoer.to_json())
    map_source = GeoJSONDataSource(geojson=shp_afwatr.to_json())
    
    # -----------------------------------------------
    # Plot map
    
    p1 = figure(title='Afvoer (laterals) T=1 zomer: RUPROF input vs. Hydromedah [m3/s]', height=350, width=820)
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_poly  = p1.patches(fill_alpha=1,
              fill_color={'field': metric, 'transform': color},
              line_color='black', line_width=0.5, source=map_source)
    color_bar = ColorBar(color_mapper=color,title=metric)
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@CODE'),
                ('Label', '@'+metric)]
    hover = HoverTool(renderers=[map_poly], tooltips=tooltips) 
    p1.add_tools(hover)  
    p1.patches(fill_alpha=0,fill_color='white',line_color='black', line_width=1.5, source=geo_source)
    
    
    output_file(paths.fightml)
    save(p1)
    
def fun_T1(x=None):
    # waarde voor T=1 afleiden
    x = x[np.isnan(x)==0]
    x = np.sort(x)[::-1]
    N = len(x)
    i = np.arange(1,N+1)
    p = i/(N+1)
    T = 1/p
    
    if N > 0:
        interp_func = interp1d(T,x)
        x_T1        = interp_func(365)
    else:
        x_T1 = -999
    
    return np.round(float(x_T1),2)
    
def add_statistic(output, data, threshold=None, percentile=None, name=None):
    # functie om statistieken toe te voegen aan de shapefile
    data = abs(data)
    
    statcol = data.apply(fun_T1, axis=0)
    statcol.name = name
    output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)     
    
    return(output)

def data_2_shp(shp, df):
    ID      = df['ID'].values    
    nmax    = len(shp)
    shp_dat = np.zeros((nmax, 1))-999
    for i in range(0,len(df)):
        ind = np.where(ID[i]==shp['ID'].values)
        if len(ind)==1: 
            shp_dat[ind[0]]= df['Q_Sobek'].values[i]         
        else:
            print(ID[i], "uitzondering 2, weggelaten...")
    shp['Q_Sobek']  = shp_dat
    return shp

def lat_data(paths):
    data_file = open(paths.dat_lat,'r')
    ID   = []
    area = []
    rain = []
    sep  = []
    for line in data_file:
        data = line.split()
        ID.append(str(data[2].replace("'","")))
        area.append(float(data[-2])) #m2
        rain.append(float(data[11])) #mm/s
        sep.append(float(data[15])) #mm/s
    
    rain = np.array(rain)
    sep  = np.array(sep)
    Qlat = (rain + sep)/ 1000 * area # m3/s
    df_Qlat = pd.DataFrame({'ID': ID,'Q_Sobek': np.round(Qlat,2)})
    return df_Qlat


def main(paths, months, years):
    '''HYDROMEDAH'''
    # shapefiles inlezen
    shp_line = gpd.read_file(paths.shp_RchSeg)
    shp_line.columns = [col.strip() for col in shp_line.columns]
    shp_line.index = shp_line['ID']
    shp_line.drop(['NAME', "TYPE", 'PARENTID', 'ID_FROM', 'ID_TO', 'USERID'], axis=1, inplace=True)
    shp_afvoer= gpd.read_file(paths.shp_afvoer)
    shp_afwatr=gpd.read_file(paths.shp_afwatr)
    shp_node  = gpd.read_file(paths.shp_node)
    shp_node.columns = [col.strip() for col in shp_node.columns]
    
    # Data inlezen 
    # (verwerkt in Sobek HIS bestanden)    
    cases = [1,2,3,4,5,6,7,8,9,10,11,12]#
    subhislist = []
    for case in cases:
        pad_to_his = os.path.join(paths.dat_Sobek, str(case), 'reachseg.his')
        reachseg = hkv.read_his.ReadMetadata(pad_to_his)    
        reach_segments = reachseg.DataFrame()
        reach_segments = reach_segments.iloc[32:,:]
        locs= reachseg.GetLocations()
        params= reachseg.GetParameters()
        reach_segments = reach_segments.loc[:,params[params.index('Discharge mean(m³/s)')]]
        reach_segments.head()
        subhislist.append(reach_segments)
    all_data = pd.concat(subhislist)
    all_data['MONTH'] = all_data.index.month
    all_data['YEAR'] = all_data.index.year
    
    # Statistiek berekenen
    sub_data = all_data[all_data['MONTH'].isin(months)]
    sub_data = all_data[all_data['YEAR'].isin(years)]
    shp_HYDROMEDAH = add_statistic(shp_line, sub_data, name='Q_HYDROMEDAH')
    
    # Data overzetten naar punten (spatial join)
    shp_HYDROMEDAH = shp_node.sjoin_nearest(shp_HYDROMEDAH, how="left", distance_col="Distances")
    
    # check distance
    ind = np.where(shp_HYDROMEDAH['Distances'].values>0.01)[0]
    if len(ind)>0: 
        print('Check file for nodes with distance > 0.01 m! : ')
        print(paths.shp_HYDROMEDAH)
    shp_HYDROMEDAH.to_file(paths.shp_HYDROMEDAH)
    
    '''RUPROF'''    
    # RUPROF lateral data inlezen
    df_Qlat = lat_data(paths)
    
    # RUPROF lateral data koppelen aan shapefile
    shp_RUPROF = data_2_shp(shp_node, df_Qlat)
    
    '''Verschil berekenen'''
    shp_comp = shp_RUPROF.copy()    
    shp_comp = shp_comp.iloc[np.where(shp_comp['Q_Sobek']!=-999)[0]]
    dat_HYDROMEDAH = np.zeros(len(shp_comp))-999
    for q in range(0,len(shp_comp)):
        ind = np.where(shp_comp['ID'].values[q]==shp_HYDROMEDAH['ID_left'].values)[0]
        dat_HYDROMEDAH[q] = shp_HYDROMEDAH['Q_HYDROMEDAH'].values[ind]
    shp_comp['Q_HYDROMEDAH'] = dat_HYDROMEDAH
    shp_comp['dQ'] = shp_comp['Q_Sobek'] - shp_comp['Q_HYDROMEDAH']
    
    # save file    
    shp_comp.rename(columns ={'ID':'lat_ID'}).to_file(paths.shp_laterals)    
    shp_RUPROF.to_file(paths.shp_RUPROF)   
    
    # add values to shapefile
    shp_comp['CODE'] = [i.replace('drain','') for i in shp_comp.ID.values]
    shp_afwatr  = shp_afwatr.merge(shp_comp[['CODE','dQ']], on='CODE', how='left')
    
    # Data plotten
    png(paths,shp_afwatr,shp_afvoer,shp_comp)
    bokeh(paths, shp_afvoer, shp_afwatr)