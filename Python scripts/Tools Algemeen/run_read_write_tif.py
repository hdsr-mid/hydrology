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
    
    return arr, ds

def write_geotiff(filename, arr, in_ds):
    if arr.dtype == np.float32:
        arr_type = gdal.GDT_Float32
    else:
        arr_type = gdal.GDT_Int32

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(filename, arr.shape[1], arr.shape[0], 1, arr_type)
    out_ds.SetProjection(in_ds.GetProjection())
    out_ds.SetGeoTransform(in_ds.GetGeoTransform())
    band = out_ds.GetRasterBand(1)
    band.WriteArray(arr)
    band.FlushCache()
    band.ComputeStatistics(False)

if __name__ == "__main__":
    file = 'file.tif'
    data1_arr, data1_ds = read_geotiff(file)
    data2_arr, data2_ds = read_geotiff(file)
    
    data3 = np.where(data1_arr != data2_arr, 1, 0)
    
    write_geotiff("file_changed.tif", data3, data1_ds)
    
    plt.subplot(311)
    plt.imshow(data1_arr)
    
    plt.subplot(312)
    plt.imshow(data2_arr)
    
    plt.subplot(313)
    plt.imshow(data3)
    
    plt.savefig('fig.png')
