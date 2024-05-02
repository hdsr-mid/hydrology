# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 16:20:18 2024

@author: petra
"""


import numpy as np
from bokeh.plotting import figure, save
from bokeh.models import CustomJS, ColumnDataSource, TapTool
from bokeh.layouts import layout
from bokeh.plotting import figure, save
from bokeh.models import ColumnDataSource, GeoJSONDataSource
from bokeh.models import LinearColorMapper
from bokeh.palettes import Viridis6 as palette
import geopandas as gpd
import pandas as pd
import numpy as np
from bokeh.models import HoverTool 
import warnings
warnings.filterwarnings("ignore")

if __name__ == "__main__":
    # -----------------------------------------------
    # Get data
    points = gpd.read_file('points.gpkg', layer = 'dmnodes')
    points['population'] = np.random.rand(len(points))
    points = points.drop('geometry', axis=1).copy()
    
    gdf  = gpd.read_file('polygons.gpkg')
    geo_source = GeoJSONDataSource(geojson=gdf.to_json())
    
    # -----------------------------------------------
    # Plot map
    p_map = figure(title="My first interactive plot!", height=300)
    p_map.patches('xs', 'ys', fill_alpha=0.7,
              fill_color={'field': 'code', 'transform': LinearColorMapper(palette=palette)},
              line_color='black', line_width=0.5, source=geo_source)
    map_source = ColumnDataSource(points)
    map_points = p_map.circle('x', 'y', source=map_source, color='red', size=10)
    tooltips = [('ID', '@id'),
                ('Label', '@population')]
    hover = HoverTool(renderers=[map_points], tooltips=tooltips) 
    p_map.add_tools(hover) 
    
    # -----------------------------------------------
    # Get timeseries
    time = pd.date_range(start='1/1/2018', end='1/01/2024')
    xx = np.linspace(0,10,len(time))
    src_dict = {}
    for s in range(0,len(points)):
        ID = points.iloc[s]['id']
        src_dict[ID] = {'x':time,'y':np.sin(xx)*np.random.rand(1)}
        
    # -----------------------------------------------
    # Tapping...
    #start with initial state for columndatasource
    ID_0 = points.iloc[0]['id']
    src=ColumnDataSource(data=src_dict[ID_0])
    p_TS = figure(title=ID_0, height=150, y_range=(0,1), x_axis_type='datetime')
    line_rend = p_TS.line(x='x', y='y',line_color='blue',source=src)
    
    #create the TapTool, add it to the bar figure, and enable it on initialization
    tap = TapTool(renderers=[map_points])
    p_map.add_tools(tap)
    p_map.toolbar.active_tap = tap
    
    cb = CustomJS(args=dict(map_source=map_source,src_dict=src_dict,src=src,p_TS=p_TS,line_rend=line_rend)
                  ,code='''
                  var sel_bar_i = map_source.selected.indices[0]
                  var line_id = map_source.data['id'][sel_bar_i]
                  src.data = src_dict[line_id]
                  p_TS.title.text=line_id
                  src.change.emit()
                  p_TS.change.emit()
                  ''')
    #apply this callback to trigger whenever the indices of the selected glyph change
    map_source.selected.js_on_change('indices',cb)             
    
    save(layout([p_map,p_TS]),'fig_clickable_hover.html')
