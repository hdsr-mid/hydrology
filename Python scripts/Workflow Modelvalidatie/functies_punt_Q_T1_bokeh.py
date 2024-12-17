# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:26:30 2024

@author: PetraH
"""

import numpy as np
import bokeh.plotting 
from bokeh.models import CustomJS, ColumnDataSource, GeoJSONDataSource, TapTool, LinearColorMapper, ColorBar, MultiLine, HoverTool
from bokeh.layouts import layout
from bokeh.transform import transform
import geopandas as gpd
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
        
def get_data(df_data_WL, shp_points, ids):
    ID_data  = shp_points.iloc[ids]['ID_left']
    ID_model = shp_points.iloc[ids]['ID_right']
    
    df_data = df_data_WL[['time',ID_data]]
    
    tijd  = df_data.time.values
    OBS   = df_data[ID_data].values    
    MOD_T1= shp_points.iloc[ids]['Model_T1'] + np.zeros(len(OBS))
    OBS_T1= shp_points.iloc[ids]['Data_T1'] + np.zeros(len(OBS))
    
    tijd  = np.reshape(tijd,(len(tijd),1))
    OBS   = np.reshape(OBS,(len(OBS),1))
    data  = [tijd, MOD_T1, OBS, OBS_T1]
    return data

def fun_create_dict(shp_points, df_data):
    # Cerate timeseries dictionary
    src_dict          = {}
    var               = ['tijd','MOD_T1','OBS','OBS_T1']
    
    # Maximum number of rows with the same coordinates
    smax = np.max(shp_points.groupby(['X','Y']).size().reset_index().rename(columns={0:'count'})['count'])
    
    for point in range(0,len(shp_points)):    
        # Extract data
        FID      = shp_points.iloc[point]['FID']
                    
        # combine timeseries at the same location
        Xi      = shp_points.iloc[point]['X']
        Yi      = shp_points.iloc[point]['Y']
        idx     = np.where((shp_points['X']==Xi) & (shp_points['Y']==Yi))[0]        
        tijd    = [get_data(df_data, shp_points, idx[ids])[0] if ids < len(idx) else get_data(df_data, shp_points, idx[0])[0] for ids in range(0,smax)]     
        MOD_T1  = [get_data(df_data, shp_points, idx[ids])[1] if ids < len(idx) else get_data(df_data, shp_points, idx[0])[1] for ids in range(0,smax)]     
        OBS     = [get_data(df_data, shp_points, idx[ids])[2] if ids < len(idx) else get_data(df_data, shp_points, idx[0])[2] for ids in range(0,smax)]     
        OBS_T1  = [get_data(df_data, shp_points, idx[ids])[3] if ids < len(idx) else get_data(df_data, shp_points, idx[0])[3] for ids in range(0,smax)]     
            
        data  = [tijd, MOD_T1, OBS, OBS_T1]
                     
        # Create dictionary
        src_dict[FID] = dict(zip(var, data))
        
    
    return src_dict


def bokeh_fig(paths, shp_afvoer, shp_points, df_data, vmin, vmax):  
    metric = 'dQ'
    
    # -----------------------------------------------
    # Get data
    shp_reach  = gpd.read_file(paths.shp_reach)
    shp_reach['x'] = shp_reach.apply(getLineCoords, geom='geometry', coord_type='x', axis=1)
    shp_reach['y'] = shp_reach.apply(getLineCoords, geom='geometry', coord_type='y', axis=1)
    shp_points = shp_points[shp_points[metric].replace(np.nan,-999)!=-999]
    df_data    = df_data.replace(-999,np.nan)
    geo_source = GeoJSONDataSource(geojson=shp_afvoer.to_json())
    
    shp_points  = shp_points.drop('geometry', axis=1).copy()
    shp_points['FID'] = shp_points.reset_index().index.values.astype(str)
    
    # -----------------------------------------------
    # Plot map: with points
    metric_abs    = abs(shp_points[metric])
    metric_norm   = metric_abs/np.max(metric_abs)
    shp_points['msize'] = np.where(abs(shp_points[metric])<0.05, 8, 8+metric_norm *20)     
    shp_points    = shp_points.rename(columns ={'X_left':'X','Y_left':'Y'})
    
    p1 = bokeh.plotting.figure(title='Debiet T=1 zomer: RUPROF vs. FEWS-WIS [m3/s]', height=350, width=820)
    p1.patches(fill_alpha=1,fill_color='white',line_color='black', line_width=0.5, source=geo_source)
    p1.multi_line('x', 'y', source=ColumnDataSource(shp_reach.drop('geometry', axis=1)), color='blue', line_width=0.2)
    map_source = ColumnDataSource(shp_points)    
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_points = p1.scatter('X', 'Y', source=map_source,color=transform(metric, color), size='msize')
    color_bar = ColorBar(color_mapper=color,title='dQ = model (RUPROF) - data (FEWS-WIS) [m3/s]')
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@ID_left'),
                ('Label', '@'+metric+'{0.0}'+'m3/s'),
                ('OBS', '@Data_T1{0.0}'+'m3/s'),
                ('MOD', '@Model_T1{0.0}'+'m3/s')]
    hover = HoverTool(renderers=[map_points], tooltips=tooltips) 
    p1.add_tools(hover)  
    
    # -----------------------------------------------
    # Plot map: without points
    p2 = bokeh.plotting.figure(title='Op basis van maximale absolute waarde per afvoergebied...', height=350, width=700)
    map_poly = p2.patches(fill_alpha=1,
              fill_color={'field': metric, 'transform': color},
              line_color='black', line_width=0.5, source=geo_source)
    tooltips = [('Code', '@NAAM'),
                ('Label', '@'+metric+'{0.0}'+'m')]
    hover = HoverTool(renderers=[map_poly], tooltips=tooltips) 
    p2.add_tools(hover)  
    
    # -----------------------------------------------
    # Merge dictionary: combine at the same location
    src_dict = fun_create_dict(shp_points, df_data)
    
    # Tapping...
    #start with initial state for columndatasource
    ID_0  = shp_points.iloc[0]['FID']
    src   = ColumnDataSource(data=src_dict[ID_0])
    p3    = bokeh.plotting.figure(title='LOC_ID ' +  shp_points.iloc[0]['ID_left'], height=300, width=1400, x_axis_type='datetime', x_axis_label="time",y_axis_label="Q [m3/s]") # , y_range=(0,1)
    glyph = MultiLine(xs="tijd", ys="OBS", line_color="black")
    p3.add_glyph(src, glyph)
    glyph = MultiLine(xs="tijd", ys="MOD_T1", line_color="red")
    p3.add_glyph(src, glyph)
    glyph = MultiLine(xs="tijd", ys="OBS_T1", line_color="black")
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
                  p3.title.text= 'LOC_ID ' + map_source.data['ID_left'][sel_bar_i]
                  src.change.emit()
                  p3.change.emit()
                  ''')
    
    
    #apply this callback to trigger whenever the indices of the selected glyph change
    map_source.selected.js_on_change('indices',cb)   
    
    
    grid = layout([[p1, p2],[p3,],])
    bokeh.plotting.show(grid, notebook_handle=True) 
    bokeh.plotting.save(grid,paths.fightml)  
        