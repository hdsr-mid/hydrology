# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 13:42:55 2024

@author: PetraH
"""

import numpy as np
import geopandas as gpd
import warnings
warnings.filterwarnings("ignore")
import bokeh.plotting
from bokeh.models import GeoJSONDataSource, LinearColorMapper, ColorBar
from bokeh.layouts import layout
from bokeh.models import HoverTool 
bokeh.plotting.output_notebook()

def fun_centroid_points(shp):
    shp = shp.copy()
    shp['geom']      = shp['geometry']
    shp['geometry']  = shp.centroid
    shp              = shp.set_geometry('geometry')    
    return shp    

def peil(shp,cols):
    # get streefpeil
    cc = ['ZP','WP','VP','BP','OP']
    Peil = []
    for p in range(0,len(shp)):
        peil = ', '.join([cc[c] + ': ' + str(shp.iloc[p][cols[c]]) for c in range(0,len(cols))])
        Peil = np.append(Peil,peil)
    
    return Peil

def bokeh_fig(paths, shp,shp_diff):  
    metric = 'diff'
    
    # -----------------------------------------------
    # Get data
    map_source1 = GeoJSONDataSource(geojson=shp.to_json())
    map_source3 = GeoJSONDataSource(geojson=shp_diff.to_json())
    
    # -----------------------------------------------
    # Plot map    
    p1 = bokeh.plotting.figure(title='Streefpeil: BR vs. Peilbesluit [m NAP]', height=350, width=820)
    color = LinearColorMapper(palette = 'Turbo256', low = 0, high = 1)
    map_poly = p1.patches(fill_alpha=1,fill_color={'field': metric, 'transform': color},line_color='black', line_width=0.5, source=map_source1)
    color_bar = ColorBar(color_mapper=color,title='Numerieke verschillen [0: nee, 1: ja]')
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID (BR)', '@CODE'), ('ID (PS)', '@CODE_PS'),('BR', '@BR'),('PS', '@PS')]
    hover = HoverTool(renderers=[map_poly], tooltips=tooltips) 
    p1.add_tools(hover)  
    
    p2 = bokeh.plotting.figure(title='Streefpeil: BR vs. Peilbesluit', height=350, width=820)
    p2.patches(fill_alpha=1, fill_color = 'grey',line_color='black', line_width=0.5, source=map_source1, legend_label='BR')
    p2.patches(fill_alpha=0, fill_color='white',line_color='red', line_width=0.5, source=map_source3, legend_label='Peilbesluit afwijkend van BR')
    p2.legend.location = "top_right"
    
    grid = layout([[p1, p2],])
    bokeh.plotting.show(grid, notebook_handle=True)
    bokeh.plotting.save(grid,paths.fightml)
    
def main(paths):
    # Load data
    shp_streefpeil     = gpd.read_file(paths.shp_streefpeil)
    shp_peilbesluit    = gpd.read_file(paths.shp_peilbesluit)
    
    # columns
    col_streefpeil  = ['ZOMERPEIL','WINTERPEIL','VASTPEIL','FLEXIBEL_B','FLEXIBEL_O']
    col_peilbesluit = ['WS_ZP','WS_WP','WS_VP','WS_BP','WS_OP']
    
    # get peil
    shp_streefpeil['BR'] = peil(shp_streefpeil,col_streefpeil)
    shp_peilbesluit['PS']= peil(shp_peilbesluit,col_peilbesluit)
    
    shp_streefpeil  = shp_streefpeil[['CODE','BR','geometry']]
    shp_peilbesluit = shp_peilbesluit[['CODE','PS','geometry']]
    
    # numeriek verschil
    points_streefpeil     = fun_centroid_points(shp_streefpeil)
    points_peilbesluit    = fun_centroid_points(shp_peilbesluit)
    shp_joined            = gpd.sjoin_nearest(points_streefpeil, points_peilbesluit[['CODE','PS','geometry']], distance_col="distances", how="left")
    shp_joined['CODE']    = shp_joined['CODE_left']
    shp_joined['CODE_PS'] = shp_joined['CODE_right']
    shp_streefpeil        = shp_streefpeil.merge(shp_joined[['CODE','CODE_PS','PS','distances']], on='CODE', how='left')
    shp_streefpeil['diff']= np.where(shp_streefpeil['BR']==shp_streefpeil['PS'],0,1)
    
    # ruimtelijk verschil    
    difference           = shp_peilbesluit.difference(shp_streefpeil, align=True).reset_index()
    difference.columns   = ['CODE','geometry']
    difference           = difference.set_geometry('geometry')
    difference['area']   = difference.geometry.area
    difference['length'] = difference.geometry.length
    difference['width']  = difference.area/difference.length
    difference           = difference.iloc[np.where(difference.width>2)[0]]
    
    # create figure
    bokeh_fig(paths, shp_streefpeil,difference)
    
    # save files
    shp_streefpeil.to_file(paths.shp_out)
    difference.to_file(paths.shp_diff)
