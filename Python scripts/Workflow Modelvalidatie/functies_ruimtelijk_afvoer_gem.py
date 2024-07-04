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
from bokeh.models import ColumnDataSource, GeoJSONDataSource, LinearColorMapper, ColorBar, Range1d
from bokeh.layouts import layout
from bokeh.transform import transform
from bokeh.models import HoverTool 



from bokeh.models import Whisker
from bokeh.sampledata.autompg2 import autompg2
from bokeh.transform import factor_cmap


vmin   = -2
vmax   = 2

# voor de boxplot en histogram
dmin  = 0
dmax  = 25
nbins = dmax-dmin

def png(paths,shp_afvoer,shp_comp):
    # Smaller marker size for when the difference is insignificant
    
    plt.figure(figsize=(10,5))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    ax = plt.subplot2grid((2, 1), (0, 0), rowspan = 2)
    shp_comp.plot(ax=ax, column = 'dQ', vmin = vmin, vmax = vmax, cmap = 'jet', edgecolor='face', linewidth=0, alpha = 1, legend=True, legend_kwds={"label": "Verschil: Satelliet - Hydromedah [mm/d]"})
    shp_afvoer.plot(ax=ax, color= "none", facecolor = "none", edgecolor='black')
    plt.xlim([109000,170000])
    plt.ylim([438000,468000])
    plt.grid()
    
    plt.savefig(paths.fig, dpi=300)
    plt.show()
    plt.close()
    
    '''-----------------------------------------------------------------------------------------------'''
    shp_comp = shp_comp.rename(columns ={'Qsat':'Satelliet','Q_HYDROMEDAH':'HYDROMEDAH'})
                                        
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.9)
    
    
    ax = plt.subplot2grid((2, 2), (0, 0))
    shp_comp.boxplot(column=['Satelliet','HYDROMEDAH'])
    plt.ylim([dmin,dmax])
    plt.ylabel('Q [mm/d]')
    
    ax = plt.subplot2grid((2, 2), (1, 0))
    shp_comp['Satelliet'].hist(bins=nbins, range=(dmin,dmax))
    plt.title('Satelliet data')
    plt.xlabel('Q [mm/d]')
    plt.xlim([dmin,dmax])
    plt.ylabel('Nr bins')
    
    ax = plt.subplot2grid((2, 2), (1, 1))
    shp_comp['HYDROMEDAH'].hist(bins=nbins, range=(dmin,dmax))
    plt.title('HYDROMEDAH')
    plt.xlabel('Q [mm/d]')
    plt.xlim([dmin,dmax])
    plt.ylabel('Nr bins')
    
    plt.savefig(paths.fig.replace('.png','_2.png'), dpi=300)
    plt.show()
    plt.close()

def bokeh_boxplot(shp_comp):
    
    df1 = pd.DataFrame({'kind': ['Satelliet']*len(shp_comp),'value': shp_comp['Qsat']})
    df2 = pd.DataFrame({'kind': ['HYDROMEDAH']*len(shp_comp),'value': shp_comp['Q_HYDROMEDAH']})
    df  = pd.concat([df1, df2],axis=0, join="inner",ignore_index=False)     
    kinds = df.kind.unique()
    
    # compute quantiles
    qs = df.groupby("kind").value.quantile([0.25, 0.5, 0.75])
    qs = qs.unstack().reset_index()
    qs.columns = ["kind", "q1", "q2", "q3"]
    df = pd.merge(df, qs, on="kind", how="left")

    # compute IQR outlier bounds
    iqr = df.q3 - df.q1
    df["upper"] = df.q3 + 1.5*iqr
    df["lower"] = df.q1 - 1.5*iqr

    source = ColumnDataSource(df)

    p = figure(x_range=kinds, tools="", toolbar_location=None,y_axis_label="Q [mm/d]")

    # outlier range
    whisker = Whisker(base="kind", upper="upper", lower="lower", source=source)
    whisker.upper_head.size = whisker.lower_head.size = 20
    p.add_layout(whisker)

    # quantile boxes
    p.vbar("kind", 0.7, "q2", "q3", source=source, line_color="black")
    p.vbar("kind", 0.7, "q1", "q2", source=source, line_color="black")

    # outliers
    outliers = df[~df.value.between(df.lower, df.upper)]
    p.scatter("kind", "value", source=outliers, size=6, color="black", alpha=1)

    p.xgrid.grid_line_color = None
    p.axis.major_label_text_font_size="14px"
    p.axis.axis_label_text_font_size="12px"
    p.y_range = Range1d(dmin,dmax) 

    return p
    
def bokeh_hist(x,title):
    
    p = figure(width=670, height=400, toolbar_location=None, title=title)

    # Histogram
    bins = np.linspace(dmin, dmax, nbins)
    
    hist, edges = np.histogram(x, density=True, bins=bins)
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],fill_color="blue", line_color="white")

    p.y_range.start = 0
    p.xaxis.axis_label = "Q [mm/d]"
    p.yaxis.axis_label = "PDF(x)"

    return p

def bokeh(paths, shp_afvoer, shp_comp):  
    
    # -----------------------------------------------
    # Get data
    geo_source = GeoJSONDataSource(geojson=shp_afvoer.to_json())
    map_source = GeoJSONDataSource(geojson=shp_comp.to_json())
    
    # -----------------------------------------------
    # Plot map 
    metric = 'dQ'
    p1 = figure(title='Afvoer (laterals) gemiddeld: Satelliet vs. HYDROMEDAH [mm/d]', height=350, width=820)
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_poly  = p1.patches(fill_alpha=1,
              fill_color={'field': metric, 'transform': color},
              line_color='black', line_width=0.5, source=map_source)
    color_bar = ColorBar(color_mapper=color)
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@CODE'),
                ('Label', '@'+metric)]
    hover = HoverTool(renderers=[map_poly], tooltips=tooltips) 
    p1.add_tools(hover)  
    
    # boxplot
    p2 = bokeh_boxplot(shp_comp)
    
    # histogram
    p3 = bokeh_hist(shp_comp['Qsat'],'Satelliet data')
    p4 = bokeh_hist(shp_comp['Q_HYDROMEDAH'],'HYDROMEDAH')
    
    grid = layout([[p1,],[p2,],[p3,p4],])
    save(grid,paths.fightml)    
    
def add_statistic(output, data, threshold=None, percentile=None, name=None):
    # functie om statistieken toe te voegen aan de shapefile
    data = abs(data)
    
    statcol = data.mean(axis=0, skipna=True)
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
    ''' Satellite data'''
    # input data
    shp_afwatr = gpd.read_file(paths.shp_afwatr)
    shp_peil   = gpd.read_file(paths.shp_peil)
    df_P       = pd.read_csv(paths.df_P, delimiter=',')
    df_E       = pd.read_csv(paths.df_E, delimiter=',')
    
    # format data
    df_P   = df_P.drop(0).drop(columns='GMT+1').astype(float)
    df_E   = df_E.drop(0).drop(columns='GMT+1').astype(float)
    
    # calculate mean
    df_P  = df_P.replace(-9999,np.nan).mean(axis=0, skipna=True).reset_index()
    df_E  = df_E.replace(-9999,np.nan).mean(axis=0, skipna=True).reset_index()    
    df_P.columns = ['CODE','Psat']
    df_E.columns = ['CODE','Esat']
    
    # add values to shapefile
    shp_sat   = shp_peil[['CODE','geometry']].merge(df_P[['CODE','Psat']], on='CODE', how='left')
    shp_sat   = shp_sat.merge(df_E[['CODE','Esat']], on='CODE', how='left')
    shp_sat['Qsat'] = (shp_sat['Psat'] - shp_sat['Esat'])/365
    
    '''HYDROMEDAH'''
    # shapefiles inlezen
    shp_line = gpd.read_file(paths.shp_RchSeg)
    shp_line.columns = [col.strip() for col in shp_line.columns]
    shp_line.index = shp_line['ID']
    shp_line.drop(['NAME', "TYPE", 'PARENTID', 'ID_FROM', 'ID_TO', 'USERID'], axis=1, inplace=True)
    shp_afvoer= gpd.read_file(paths.shp_afvoer)
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
    
    # convert based on area afwateringsgebied (m3/s -> mm/d)
    shp_afwatr['area']              = shp_afwatr.area
    shp_HYDROMEDAH['CODE']          = [i.replace('drain','') for i in shp_HYDROMEDAH.ID_left.values]    
    shp_HYDROMEDAH                  = shp_HYDROMEDAH.merge(shp_afwatr[['CODE','area']], on='CODE', how='left')
    shp_HYDROMEDAH['Q_HYDROMEDAH']  = shp_HYDROMEDAH['Q_HYDROMEDAH']/shp_HYDROMEDAH['area'] * 1e3 * 3600*24
    
    
    # shapefiles opslaan
    shp_sat.to_file(paths.shp_sat)    
    shp_HYDROMEDAH.to_file(paths.shp_HYDROMEDAH)    
   
    # Data overzetten naar shp_peil (spatial join)
    shp_peil.crs       ='EPSG:28992'
    shp_HYDROMEDAH.crs ='EPSG:28992'
    shp_HYDROMEDAH     = shp_HYDROMEDAH.rename(columns ={'ID_left':'ID1','ID_right':'ID2'})[['ID1','ID2','Q_HYDROMEDAH','geometry']]
    shp_HYDROMEDAH     = gpd.sjoin(shp_HYDROMEDAH, shp_peil[['CODE','geometry']], how="left")
    shp_HYDROMEDAH     = shp_HYDROMEDAH.groupby('CODE')['Q_HYDROMEDAH'].agg(['mean']).reset_index()
    shp_HYDROMEDAH.columns = ['CODE','Q_HYDROMEDAH']
    
    # tijdelijk handmatig bestand hier gebruikt ivm verouderde peilgebieden in WIS
    shp_sat  = gpd.read_file(paths.shp_temp) 
    shp_sat['Qsat'] = shp_sat['Qsat'] / 365
    
    '''Samenvoegen & verschil berekenen'''
    shp_comp = shp_peil.merge(shp_sat.drop(columns='geometry'), on='CODE', how='left')
    shp_comp = shp_comp.merge(shp_HYDROMEDAH, on='CODE', how='left')
    shp_comp['dQ'] = np.round(shp_comp['Qsat'] - shp_comp['Q_HYDROMEDAH'],1)
    
    # print(np.min(shp_comp.dQ),np.max(shp_comp.dQ))
    # save file    
    shp_comp.to_file(paths.shp_comp)    
    
    
    # Data plotten
    png(paths,shp_afvoer,shp_comp)
    bokeh(paths, shp_afvoer,shp_comp)
    