#!/usr/bin/env python
# coding: utf-8


import os
import sys
import numpy as np
# import pandas as pd
import warnings
import geopandas as gpd
from time import sleep
warnings.filterwarnings("ignore")
from multiprocessing import Process

package_path    = os.path.abspath('../')
root            = os.getcwd()
GIS_path        = os.path.abspath(os.path.join(root,'GIS'))
temp_path       = os.path.abspath(os.path.join(root,'temp'))
sys.path.append(package_path)

from src.waterberging import WaterBerging, fun_combine_shp, fun_afronden

if __name__ == "__main__":
   
    # run per peilgebied
    WaterBerging(GIS_path=GIS_path,temp_path=temp_path, export_path=GIS_path)
    
    # combineer resultaten
    fun_combine_shp(GIS_path=GIS_path,temp_path=temp_path, export_path=GIS_path)
    
    # afronden
    fun_afronden(GIS_path=GIS_path,temp_path=temp_path, export_path=GIS_path)






