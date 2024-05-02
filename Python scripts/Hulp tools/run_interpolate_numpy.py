# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 16:05:45 2023

@author: PetraH
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

if __name__ == "__main__":  
    
    #Let's create some random  data
    data_orig = np.random.random_integers(0,10,(10,10)).astype(float)
    #values grater then 7 goes to np.nan
    data_orig[data_orig>7] = np.nan
    
    x = np.arange(0, data_orig.shape[1])
    y = np.arange(0, data_orig.shape[0])
    
    #mask invalid values
    array = np.ma.masked_invalid(data_orig)
    xx, yy = np.meshgrid(x, y)
    
    #get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                              (xx, yy),
                                 method='linear')