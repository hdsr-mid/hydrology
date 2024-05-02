# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 12:47:59 2024

@author: petra
"""

from bokeh.plotting import figure, save
from bokeh.models import ColumnDataSource, GeoJSONDataSource
from bokeh.models import LinearColorMapper
from bokeh.palettes import Viridis6 as palette

import geopandas as gpd
import numpy as np
from bokeh.models import HoverTool 
import warnings
warnings.filterwarnings("ignore")

if __name__ == "__main__":
    points = gpd.read_file('points.gpkg', layer = 'dmnodes')
    points['population'] = np.random.rand(len(points))
    points = points.drop('geometry', axis=1).copy()
    
    gdf  = gpd.read_file('polygons.gpkg')
    geo_source = GeoJSONDataSource(geojson=gdf.to_json())
    
    # Initialize the plot (p) and give it a title
    p = figure(title="My first interactive plot!")
    p.patches('xs', 'ys', fill_alpha=0.7,
              fill_color={'field': 'code', 'transform': LinearColorMapper(palette=palette)},
              line_color='black', line_width=0.5, source=geo_source)
    r = p.circle('x', 'y', source=ColumnDataSource(points), color='red', size=10)
    tooltips = [('ID', '@id'),
                ('Label', '@population')]
    hover = HoverTool(renderers=[r], tooltips=tooltips) 
    p.add_tools(hover) 
    
    # Save the plot by passing the plot -object and output path
    outfp = "fig_hover.html"
    save(obj=p, filename=outfp)