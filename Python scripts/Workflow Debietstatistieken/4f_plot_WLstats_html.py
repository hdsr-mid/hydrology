# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 13/04/2026
virtual environment used: geo-env + bokeh (apart installeren; code: conda install conda-forge::bokeh)

"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely import get_coordinates
from shapely.geometry import MultiLineString

import bokeh.plotting
from bokeh.models import CustomJS, ColumnDataSource, GeoJSONDataSource
from bokeh.models import LinearColorMapper, ColorBar, Line, TapTool
from bokeh.layouts import layout
from bokeh.transform import transform
from bokeh.models import HoverTool 
bokeh.plotting.output_notebook()

import warnings
warnings.filterwarnings("ignore")


class general():
    root = r'H:/Team_Kennis_OSA/Hydrologische_Informatieproducten/4. Debietstatistieken (30042026)'
    path_all            = root + '/03_Output/netcdf/DFM_all_Waterstand.nc'
    outdir              = root + '/03_Output'
    path_WLstats        = root + '/03_Output/gpkg/Waterstandstatistieken.gpkg'
    shp_afvoergebieden  = root + '/01_Input/Afvoergebieden_30032026.gpkg'
    shp_watergangen     = root + '/01_Input/Branches.gpkg'
   
def getLineCoords(row, geom, coord_type):
    if isinstance(row[geom], MultiLineString):
        if coord_type == 'x':
            return get_coordinates(row[geom])[:,0]
        elif coord_type == 'y':
            return get_coordinates(row[geom])[:,1]
    else:
        if coord_type == 'x':
            return list( row[geom].coords.xy[0] )
        elif coord_type == 'y':
            return list( row[geom].coords.xy[1] )
        
   
def fun_create_dict(gdf_sel,xds):
    # Create timeseries dictionary
    src_dict          = {}
    var               = ['x','y']
    var_dat           = 'mesh1d_s1' # WL
    time              = [pd.to_datetime(t) for t in xds.time.values]
    IDs               = np.array([i.strip().replace('hdsr_wa_','') for i in xds['mesh1d_node_id'].values.astype(str)])
    for iD in gdf_sel.CODE.values:    
        idx_mesh1d = int(np.where(iD==IDs)[0][0])
        data  = xds[var_dat].sel(mesh1d_nNodes=idx_mesh1d).values
        data  = np.round(data,2)
    
        data  = [time, data]
                     
        # Create dictionary
        src_dict[iD] = dict(zip(var, data))
        
    return src_dict

def fun_plot_html(ig, xds,gdf, shp_afvoer):
    # select region
    gebied = shp_afvoer.iloc[ig]
    CODE   = gebied.CODE
    NAAM   = gebied.NAAM
    gdf_sel = gdf.clip(gebied.geometry)
    if len(gdf_sel)==0: return
    
    # load new data
    watergangen = gpd.read_file(general.shp_watergangen).clip(gebied.geometry)
    watergangen['x'] = watergangen.apply(getLineCoords, geom='geometry', coord_type='x', axis=1)
    watergangen['y'] = watergangen.apply(getLineCoords, geom='geometry', coord_type='y', axis=1)
    
    # create dict
    src_dict = fun_create_dict(gdf_sel,xds)
    tmax     = len(xds.time.values)
    for iD in gdf_sel.CODE.values: src_dict[iD]['T1'] = gdf_sel.loc[gdf_sel.CODE==iD]['WL_T1'].values*np.ones(tmax,)
    for iD in gdf_sel.CODE.values: src_dict[iD]['gem'] = gdf_sel.loc[gdf_sel.CODE==iD]['WLgem'].values*np.ones(tmax,)
    for iD in gdf_sel.CODE.values: src_dict[iD]['10%'] = gdf_sel.loc[gdf_sel.CODE==iD]['WL_p10'].values*np.ones(tmax,)
    for iD in gdf_sel.CODE.values: src_dict[iD]['90%'] = gdf_sel.loc[gdf_sel.CODE==iD]['WL_p90'].values*np.ones(tmax,)
    
    # plot WL
    geo_source  = GeoJSONDataSource(geojson=shp_afvoer[ig:ig+1].to_json())
    map_source  = ColumnDataSource(gdf_sel.drop('geometry', axis=1))
    map_source2 = ColumnDataSource(watergangen.drop('geometry', axis=1))
    color       = LinearColorMapper(palette = 'Turbo256', low = -5, high = 5)
    
    p1 = bokeh.plotting.figure(title='Debiet', height=350, width=820, match_aspect=True)
    p1.patches(fill_alpha=0.7,fill_color='grey',line_color='black', line_width=0.5, source=geo_source)
    map_line = p1.scatter('X', 'Y', source=map_source, color=transform('WL_T1', color))
    p1.multi_line('x', 'y', source=map_source2, color='grey')
    color_bar = ColorBar(color_mapper=color,title='WL [m]')
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@CODE'),
                ('T1', '@WL_T1{0.01}'+'m'),
                ('P10', '@WL_p10{0.01}'+'m'),
                ('P90', '@WL_p90{0.01}'+'m')]
    p1.add_tools(HoverTool(renderers=[map_line], tooltips=tooltips))  
    
    # Tapping...
    #start with initial state for columndatasource
    ID_0  = list(src_dict.keys())[0]
    src   = ColumnDataSource(data=src_dict[ID_0])
    p3    = bokeh.plotting.figure(title=ID_0, height=300, width=1400, x_axis_type='datetime', x_axis_label="time",y_axis_label="WL [m]") # , y_range=(0,1)
    glyph = Line(x="x", y='y', line_color="black")
    p3.add_glyph(src, glyph)
    p3.add_glyph(src, Line(x="x", y='T1', line_color="red"))
    p3.add_glyph(src, Line(x="x", y='10%', line_color="blue", line_dash = 'dotted'))
    p3.add_glyph(src, Line(x="x", y='90%', line_color="blue", line_dash = 'dotted'))
    p3.line(x=src_dict[ID_0]['x'][0],y=src_dict[ID_0]['y'][0], line_color="black",legend_label ='WL model')
    p3.line(x=src_dict[ID_0]['x'][0],y=src_dict[ID_0]['T1'][0], line_color="red",legend_label ='T1')
    p3.line(x=src_dict[ID_0]['x'][0],y=src_dict[ID_0]['10%'][0], line_color="blue",legend_label ='10%')
    p3.line(x=src_dict[ID_0]['x'][0],y=src_dict[ID_0]['90%'][0], line_color="blue",legend_label ='90%')
    p3.legend.location = "top_left"
    
    #create the TapTool, add it to the bar figure, and enable it on initialization
    tap = TapTool(renderers=[map_line])
    p1.add_tools(tap)
    p1.toolbar.active_tap = tap
    
    cb = CustomJS(args=dict(map_source=map_source,src_dict=src_dict,src=src,p3=p3,glyph=glyph)
                  ,code='''
                  var sel_bar_i = map_source.selected.indices[0]
                  var area_id = map_source.data['CODE'][sel_bar_i]
                  p3.title.text= area_id
                  p3.change.emit()
                  src.data = src_dict[area_id]
                  src.change.emit()
                  ''')
    
    #apply this callback to trigger whenever the indices of the selected glyph change
    map_source.selected.js_on_change('indices',cb)             
    
    grid = layout([[p1],[p3,],])
    bokeh.plotting.show(grid, notebook_handle=True)
    bokeh.plotting.save(grid,os.path.join(general.outdir,'html','Waterstand',str(ig) + '_fig_WL_' + CODE + '_' + NAAM+'.html'))
    
    del src   


if __name__ == "__main__":
    
    # # load data
    print('Loading data...')
    gdf        = gpd.read_file(general.path_WLstats)
    xds_all    = xr.open_dataset(general.path_all)
    shp_afvoer = gpd.read_file(general.shp_afvoergebieden)
    
    # plot data
    print('Plotting...')
    for ig in range(0,len(shp_afvoer)):
        print(ig,len(shp_afvoer))
        fun_plot_html(ig, xds_all,gdf, shp_afvoer)
    
    xds_all.close()
    print('Done!')
    
    