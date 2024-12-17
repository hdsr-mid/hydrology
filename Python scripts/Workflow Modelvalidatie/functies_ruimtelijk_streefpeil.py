# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 13:47:02 2024

@author: PetraH
"""


import numpy as np
import pandas as pd
import geopandas as gpd
import warnings
warnings.filterwarnings("ignore")

import bokeh.plotting
from bokeh.models import ColumnDataSource, GeoJSONDataSource, LinearColorMapper, ColorBar
from bokeh.layouts import layout
from bokeh.transform import transform
from bokeh.models import HoverTool 
bokeh.plotting.output_notebook()

vmin = -1
vmax = 1

def bokeh_fig(paths, shp_afvoer, shp_points):  
    metric = 'dWL'
    
    # -----------------------------------------------
    # Get data
    shp_points = shp_points[shp_points[metric].replace(np.nan,-999)!=-999]
    shp_points = shp_points.drop('geometry', axis=1).copy()
    geo_source = GeoJSONDataSource(geojson=shp_afvoer.to_json())
    
    # -----------------------------------------------
    # Plot map: with points
    metric_abs    = abs(shp_points[metric])
    metric_norm   = metric_abs/np.max(metric_abs)
    shp_points['msize'] = np.where(abs(shp_points[metric])<0.05, 0, 1+metric_norm *20)     
    
    p1 = bokeh.plotting.figure(title='Streefpeil zomer: RUPROF (t=0) vs. BR [m NAP]', height=350, width=820)
    p1.patches(fill_alpha=1,fill_color='white',line_color='black', line_width=0.5, source=geo_source)
    map_source = ColumnDataSource(shp_points)    
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_points = p1.scatter('X', 'Y', source=map_source,color=transform(metric, color), size='msize')
    p1.scatter('X', 'Y', source=map_source,color='grey', size=0.1)
    color_bar = ColorBar(color_mapper=color,title='dWL = model (RUPROF(t=0)) - data (BR)  [m]')
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@ID'),
                ('Value', '@'+metric+'{0.0}'+'m'),
                ('BR', '@Peil{0.0}'+'m'),
                ('MOD', '@WL{0.0}'+'m')]
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
    
    grid = layout([[p1, p2],])
    bokeh.plotting.show(grid, notebook_handle=True)
    bokeh.plotting.save(grid,paths.fightml)
    
def main(paths):
    # Load data
    df_WL     = pd.read_csv(paths.df_WL_csv,skiprows=4, delimiter = ',').reset_index()
    shp_afvoer= gpd.read_file(paths.shp_afvoer)
    shp_reach = gpd.read_file(paths.shp_reach)
    shp_node  = gpd.read_file(paths.shp_node)
    df_WL.columns = ['index','ID','WL']    
    
    columns_shp_reach = [c.replace('  ','').replace('  ','').replace(' ','') for c in shp_reach.columns]
    columns_shp_node  = [c.replace('  ','').replace('  ','').replace(' ','') for c in shp_node.columns]
    shp_reach.columns = columns_shp_reach
    shp_node.columns  = columns_shp_node
    
    # Streefpeil
    shp_node=shp_node[~shp_node["ID"].str.contains('stor', na=False)] # remove stor-nodes
    Peil = np.zeros((len(shp_node),))
    for p in range(0,len(Peil)):
        label = str(shp_node['WS_LABEL'].values[p])
        if ('zp' in label):
            Peil[p] = shp_node['ZOMERPEIL'].values[p]
        elif ('vastpeil' in label):
            Peil[p] = shp_node['VASTPEIL'].values[p]
        elif ('flexibel' in label) or ('bovenpeil' in label):
            Peil[p] = shp_node['FLEXIBEL_B'].values[p] ##### check this!!        
        else:
            Peil[p] = -999            
    shp_node['Peil'] = Peil
    
    # Link Sobek data to shapefile: Waterniveau
    nmax         = len(shp_node['ID'].values)
    shp_node_WL  = np.zeros((nmax, 1))-999
    for i in range(0,len(df_WL)):
        ind = np.where((df_WL['ID'][i]==shp_node['ID'].values))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_node_WL[ind]  = df_WL['WL'].values[i]
    shp_node['WL']  = np.round(shp_node_WL,1)
    shp_node        = shp_node[shp_node['WL']!=-999]
    shp_node['dWL'] = np.round(np.where(shp_node['WL']==-999, -999  , shp_node['WL'] - shp_node['Peil']  ),1)
    shp_node['dWL'] = np.where(shp_node['Peil']==-999, -999  , shp_node['dWL'])           
    shp_node        = shp_node[shp_node['dWL']!=-999]
    shp_node.to_file(paths.shp_node_dWL)  
    
    # spatial join (nodes & shapefile)
    shp_node_joined = gpd.sjoin(shp_node[['dWL','geometry']], shp_afvoer[['CODE','geometry']], how="left")
    shp_node_joined['dWLabs'] = abs(shp_node_joined['dWL'])
    shp_node_joined = shp_node_joined.loc[shp_node_joined.groupby('CODE')['dWLabs'].idxmax().values]
    shp_afvoer      = shp_afvoer.merge(shp_node_joined[['CODE','dWL']], on='CODE', how='left')
    
    # Plot
    bokeh_fig(paths, shp_afvoer, shp_node)
    