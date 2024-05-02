# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:20:41 2023

@author: PetraH
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib.dates as mdates
from datetime import datetime
from functions import main
import warnings
warnings.filterwarnings("ignore")
import config
import functies_AHN

    
if __name__ == "__main__":
    print('Load data')
    suffices = ['init','2u_GLG','24u_GLG','96u_GLG','2u_GHG','24u_GHG','96u_GHG']    
        
    # Prepare AHN for inundation maps once
    if not os.path.exists(config.raster_clip):
        # Prepare inundationmap: Clip AHN raster to target region    
        print('Clip raster')
        functies_AHN.fun_clip_raster()
        
        # Prepare inundationmap: AHN gapfilling & set height water level (i.e., op de watergangen is de hoogte gelijk aan het waterpeil)
        print('Preprocess raster')
        functies_AHN.fun_preprocess_raster()
    
    for suffix in suffices:
        print(suffix)
        
        # Input data
        df_WL     = pd.read_csv(config.df_WL.replace('.txt','_'+suffix+'.txt'))
        
        # Time step with largest waterlevel increase
        dWL       = df_WL[df_WL.columns[2:]]
        dWL       = dWL - dWL.iloc[0]
        dWLmax    = dWL.max()
        t_WLmax   = int(np.median(np.argmax(dWL.values,axis=0)))
        
        if 'init' in suffix:
            # initial conditions
            main(0, suffix)
        else:
            # per hour    
            # main(0, suffix)
            tmax = len(df_WL['Location:'].values)
            for t in tqdm(range(0,tmax)):
                main(t, suffix)