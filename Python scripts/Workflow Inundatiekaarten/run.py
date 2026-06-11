# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:20:41 2023

@author: PetraH
"""

import functions
import functies_AHN
import warnings
warnings.filterwarnings("ignore")


if __name__ == "__main__":
    print('Load data')
    
    # Clip AHN raster to target region    
    print('Clip raster')
    functies_AHN.fun_clip_raster()
    
    # AHN gapfilling & set height water level (i.e., op de watergangen is de hoogte gelijk aan het waterpeil)
    print('Preprocess raster')
    functies_AHN.fun_preprocess_raster()
        
    # Estimate inundation depth
    print('Estimate inundation depth')
    functions.main()
    
    