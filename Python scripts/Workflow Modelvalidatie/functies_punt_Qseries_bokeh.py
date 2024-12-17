# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:26:33 2024

@author: PetraH
"""

import numpy as np
import bokeh.plotting
from bokeh.models import CustomJS, ColumnDataSource, GeoJSONDataSource, TapTool, LinearColorMapper, ColorBar, MultiLine
from bokeh.layouts import layout
from bokeh.transform import transform
import geopandas as gpd
from bokeh.models import HoverTool 
from shapely.geometry import MultiLineString
import warnings
warnings.filterwarnings("ignore")
bokeh.plotting.output_notebook()

def getLineCoords(row, geom, coord_type):
    if isinstance(row[geom], MultiLineString):
        empty_l = []
        return empty_l
    else:
        if coord_type == 'x':
            return list( row[geom].coords.xy[0] )
        elif coord_type == 'y':
            return list( row[geom].coords.xy[1] )
        
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
    shp_reach  = gpd.read_file(paths.shp_reach)
    shp_reach['x'] = shp_reach.apply(getLineCoords, geom='geometry', coord_type='x', axis=1)
    shp_reach['y'] = shp_reach.apply(getLineCoords, geom='geometry', coord_type='y', axis=1)
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
        shp_points['msize'] = np.where(abs(shp_points[metric])<1, 8, 8+metric_norm *20)     
    else:
        shp_points['msize'] = 10
    
    p1 = bokeh.plotting.figure(title='Debiet tijd series: Sobek vs. FEWS-WIS [%]', height=350, width=820)
    p1.patches(fill_alpha=1,fill_color='white',line_color='black', line_width=0.5, source=geo_source)
    p1.multi_line('x', 'y', source=ColumnDataSource(shp_reach.drop('geometry', axis=1)), color='blue', line_width=0.2)
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_source1 = ColumnDataSource(shp_points)
    map_source2 = ColumnDataSource(shp_points[shp_points[metric]!=-999])
    map_points1 = p1.scatter('X', 'Y', source=map_source1,color='grey', size=5)
    map_points2 = p1.scatter('X', 'Y', source=map_source2,color=transform(metric, color), size='msize')
    color_bar = ColorBar(color_mapper=color,title='Qmod/Qobs [%]')
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@ID_left' + '/' + '@IRIS_ID'),
                ('Value', '@'+metric+'{0.}'+'%'),
                ('OBS', '@Q_OBS{0.00}'+'m3/s'),
                ('MOD', '@Q_MOD{0.00}'+'m3/s')]
    hover = HoverTool(renderers=[map_points1], tooltips=tooltips) 
    p1.add_tools(hover)  
    
    # -----------------------------------------------
    # Plot map: without points
    p2 = bokeh.plotting.figure(title='Op basis van maximale absolute waarde per afvoergebied...', height=350, width=700)
    map_poly = p2.patches(fill_alpha=1,
              fill_color={'field': metric, 'transform': color},
              line_color='black', line_width=0.5, source=geo_source)
    tooltips = [('Code', '@NAAM'),
                ('Label', '@'+metric+'{0.}'+'%')]
    hover = HoverTool(renderers=[map_poly], tooltips=tooltips) 
    p2.add_tools(hover)  
        
    # -----------------------------------------------
    # Tapping...
    #start with initial state for columndatasource
    ID_0  = shp_points.iloc[0]['FID']
    src   = ColumnDataSource(data=src_dict[ID_0])
    p3    = bokeh.plotting.figure(title='Gemaal/ADCP ' +  shp_points.iloc[0]['IRIS_ID'], height=300, width=1400, x_axis_type='datetime', x_axis_label="time",y_axis_label="Q [m3/s]") # , y_range=(0,1)
    glyph = MultiLine(xs="t_MOD", ys="MOD", line_color="red")
    p3.add_glyph(src, glyph)
    glyph = MultiLine(xs="t_OBS", ys="OBS", line_color="black")
    p3.add_glyph(src, glyph)
    
    #create the TapTool, add it to the bar figure, and enable it on initialization
    tap = TapTool(renderers=[map_points1])
    p1.add_tools(tap)
    p1.toolbar.active_tap = tap
    
    cb = CustomJS(args=dict(map_source=map_source1,src_dict=src_dict,src=src,p3=p3,glyph=glyph)
                  ,code='''
                  var sel_bar_i = map_source.selected.indices[0]
                  var line_id = map_source.data['FID'][sel_bar_i]
                  src.data = src_dict[line_id]
                  p3.title.text= 'Gemaal/ADCP ' + map_source.data['IRIS_ID'][sel_bar_i]
                  src.change.emit()
                  p3.change.emit()
                  ''')
    
    
    #apply this callback to trigger whenever the indices of the selected glyph change
    map_source1.selected.js_on_change('indices',cb)             
    
    grid = layout([[p1, p2],[p3,],])
    bokeh.plotting.show(grid, notebook_handle=True)
    bokeh.plotting.save(grid,paths.fightml.replace('.html','_' + metric + '.html'))
    