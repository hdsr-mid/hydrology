# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:56:53 2023

@author: PetraH
"""

import os
import numpy as np
from PIL import Image
from osgeo import gdal, osr
import matplotlib.pyplot as plt
from scipy import interpolate

def read_geotiff(filename):
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    
    width  = ds.RasterXSize
    height = ds.RasterYSize
    gt     = ds.GetGeoTransform()
    x_start= gt[0] # left
    dx     = gt[1]
    x_end  = gt[0] + width*dx + height*gt[2] # right    
    y_start= gt[3] # top
    dy     = gt[5]
    y_end  = gt[3] + width*gt[4] + height*dy # bottom
    
    
    
    '''
    GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
    GT(1) w-e pixel resolution / pixel width.
    GT(2) row rotation (typically zero).
    GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
    GT(4) column rotation (typically zero).
    GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).    
    '''
    
    x = np.arange(x_start,x_end,dx)
    y = np.arange(y_start,y_end,dy)
    
    if (len(x) != width) or (len(y)!=height):
        print('Error: check dimensions')
    else:
        print('Dimensions ok')
    
    return arr, ds, x, y

if __name__ == "__main__":
    file = 'file.tif'
    data_arr, data_ds, x, y = read_geotiff(file)
    
    
    print(len(x),len(y),data_arr.shape)


