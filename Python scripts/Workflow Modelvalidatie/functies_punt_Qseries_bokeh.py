# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:26:33 2024

@author: PetraH
"""

import numpy as np
from bokeh.plotting import figure, save
from bokeh.models import CustomJS, ColumnDataSource, GeoJSONDataSource, TapTool, LinearColorMapper, ColorBar, Select, MultiLine
from bokeh.layouts import layout, gridplot, column, row
from bokeh.plotting import figure, save
from bokeh.palettes import Viridis6 as palette
from bokeh.transform import transform
import geopandas as gpd
import pandas as pd
import numpy as np
from bokeh.models import HoverTool 
import warnings
warnings.filterwarnings("ignore")

def get_data(sub_data, df_data_Q, shp_points, ids):
    ID_data  = shp_points.iloc[ids]['ID_left']
    ID_model = shp_points.iloc[ids]['ID_right']
    
    df_data = df_data_Q[['time',ID_data]]
    df_model= sub_data[ID_model].reset_index()
    df_model.columns =['time', ID_model]
    df_model = df_model.set_index('time').resample('D').mean().reset_index()
    df_data  = df_data.set_index('time').resample('D').mean().reset_index()
    
    t_MOD = df_model.time.values
    t_OBS = df_data.time.values
    MOD   = df_model[ID_model].values
    OBS   = df_data[ID_data].values    
    t_MOD = np.reshape(t_MOD,(len(t_MOD),1))
    t_OBS = np.reshape(t_OBS,(len(t_OBS),1))
    MOD   = np.reshape(MOD,(len(MOD),1))
    OBS   = np.reshape(OBS,(len(OBS),1))
    data  = [t_MOD, MOD, t_OBS, OBS]
    return data

def fun_create_dict(shp_points, sub_data, df_data_Q):
    # Cerate timeseries dictionary
    src_dict          = {}
    var               = ['t_MOD','MOD','t_OBS','OBS']
    
    # Maximum number of rows with the same coordinates
    smax = np.max(shp_points.groupby(['X','Y']).size().reset_index().rename(columns={0:'count'})['count'])
    
    for point in range(0,len(shp_points)):    
        # Extract data
        FID      = shp_points.iloc[point]['FID']
    
        # ID_data  = shp_points.iloc[point]['ID_left']
        # ID_model = shp_points.iloc[point]['ID_right']
        # data = get_data(sub_data, df_data_Q, ID_data, ID_model)
                
        # combine timeseries at the same location
        Xi    = shp_points.iloc[point]['X']
        Yi    = shp_points.iloc[point]['Y']
        idx   = np.where((shp_points['X']==Xi) & (shp_points['Y']==Yi))[0]        
        t_MOD = [get_data(sub_data, df_data_Q, shp_points, idx[ids])[0] if ids < len(idx) else get_data(sub_data, df_data_Q, shp_points, idx[0])[0] for ids in range(0,smax)]     
        MOD   = [get_data(sub_data, df_data_Q, shp_points, idx[ids])[1] if ids < len(idx) else get_data(sub_data, df_data_Q, shp_points, idx[0])[1] for ids in range(0,smax)]     
        t_OBS = [get_data(sub_data, df_data_Q, shp_points, idx[ids])[2] if ids < len(idx) else get_data(sub_data, df_data_Q, shp_points, idx[0])[2] for ids in range(0,smax)]     
        OBS   = [get_data(sub_data, df_data_Q, shp_points, idx[ids])[3] if ids < len(idx) else get_data(sub_data, df_data_Q, shp_points, idx[0])[3] for ids in range(0,smax)]     
            
        data  = [t_MOD, MOD, t_OBS, OBS]
                     
        # Create dictionary
        src_dict[FID] = dict(zip(var, data))
        
    
    return src_dict

def main(paths, shp_afvoer, shp_points, sub_data, df_data_Q, metric, vmin, vmax):    
    # -----------------------------------------------
    # Get data
    shp_points = shp_points[shp_points[metric].replace(np.nan,-999)!=-999]
    shp_points = shp_points.drop('geometry', axis=1).copy()
    shp_points['FID'] = shp_points.reset_index().index.values.astype(str)
    geo_source = GeoJSONDataSource(geojson=shp_afvoer.to_json())
    
    # Merge dictionary: combine at the same location
    src_dict = fun_create_dict(shp_points, sub_data, df_data_Q)
    
    # -----------------------------------------------
    # Plot map: with points
    metric_abs    = abs(shp_points[metric])
    metric_norm   = metric_abs/np.max(metric_abs)
    if metric == 'dQ': 
        shp_points['msize'] = np.where(abs(shp_points[metric])<0.01, 5, 5+metric_norm *20)     
    else:
        shp_points['msize'] = 10
    
    p1 = figure(title='Sobek vs. FEWS_WIS debiet data [m3/s]: ' + metric, height=350, width=820)
    p1.patches(fill_alpha=0.7,fill_color='white',line_color='black', line_width=0.5, source=geo_source)
    map_source = ColumnDataSource(shp_points)    
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_points = p1.scatter('X', 'Y', source=map_source,color=transform(metric, color), size='msize')
    color_bar = ColorBar(color_mapper=color,title=metric)
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@IRIS_ID'),
                ('Value', '@'+metric)]
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
        
    # -----------------------------------------------
    # Tapping...
    #start with initial state for columndatasource
    ID_0  = shp_points.iloc[0]['FID']
    src   = ColumnDataSource(data=src_dict[ID_0])
    p3    = figure(title='Gemaal ' +  shp_points.iloc[0]['IRIS_ID'], height=300, width=1400, x_axis_type='datetime', x_axis_label="time",y_axis_label="Q [m3/s]") # , y_range=(0,1)
    glyph = MultiLine(xs="t_MOD", ys="MOD", line_color="red")
    p3.add_glyph(src, glyph)
    glyph = MultiLine(xs="t_OBS", ys="OBS", line_color="black")
    p3.add_glyph(src, glyph)
    
    #create the TapTool, add it to the bar figure, and enable it on initialization
    tap = TapTool(renderers=[map_points])
    p1.add_tools(tap)
    p1.toolbar.active_tap = tap
    
    cb = CustomJS(args=dict(map_source=map_source,src_dict=src_dict,src=src,p3=p3,glyph=glyph)
                  ,code='''
                  var sel_bar_i = map_source.selected.indices[0]
                  var line_id = map_source.data['FID'][sel_bar_i]
                  src.data = src_dict[line_id]
                  p3.title.text= 'Gemaal ' + map_source.data['IRIS_ID'][sel_bar_i]
                  src.change.emit()
                  p3.change.emit()
                  ''')
    
    
    #apply this callback to trigger whenever the indices of the selected glyph change
    map_source.selected.js_on_change('indices',cb)             
    
    grid = layout([[p1, p2],[p3,],])
    save(grid,paths.fightml.replace('.html','_' + metric + '.html'))
    