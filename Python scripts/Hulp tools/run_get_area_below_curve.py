# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:05:16 2024

@author: PetraH
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from shapely.geometry import Polygon
import geopandas as gpd

if __name__ == "__main__":
    #create range of x-values from -4 to 4 in increments of .001
    x = np.arange(-4, 4, 0.001)
    y = -norm.pdf(x,0,1)
    x = x - x.min()
    x = np.append(x,x[0])
    y = np.append(y,y[0])
    
    # create polygon
    polygon_geom = Polygon(zip(x, y))
    polygon = gpd.GeoDataFrame(index=[0], crs='epsg:28992', geometry=[polygon_geom])   
    print(polygon.area)
