# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 13:17:46 2024

@author: PetraH
"""

import numpy as np
import os
import geopandas as gpd
import pandas as pd
from osgeo import gdal 
from shapely.geometry import Point
from shapely.geometry import LineString 
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")

def fun_filter(paths):
    
    # Input data
    gdf    = gpd.read_file(paths.shp_points_AHN_all)
    
    # Group points based on the column 'profielmet'
    grouped = gdf.groupby('profielmet')
    
    # Loop through each profielmet-code
    for name, group in tqdm(grouped):         
        '''Apply a filter based on the data availability 
        If there is no data in AHN, then the value is  around -1e22'''        
        data_min       = np.min(group[(group['punttype'] == 'A24') | (group['punttype'] == 'A28')]['profiel'].values)
        
        # Apply filter based on height difference
        if (data_min < -1e10):
            code = 0 # remove profile            
        else:
            code = 1 # keep profile
            
        # add info to gdf
        gdf.loc[(gdf['profielmet'] == name), 'filter_gaps'] = code
    
    
        '''Apply a filter based on the heights at the profile edge (insteek) and waterline: 
            height at the profile edge cannot be lower than at the waterline'''    
        # Right side
        p_insteek    = group[group['punttype'] == 'A28']
        p_water      = group[group['punttype'] == 'A22R']
        diff_R       = p_insteek['profiel'].values - p_water['profiel'].values
        
        # Left side
        p_insteek    = group[group['punttype'] == 'A24']
        p_water      = group[group['punttype'] == 'A22L']
        diff_L       = p_insteek['profiel'].values - p_water['profiel'].values
        
        
        # Apply filter based on height difference
        if (diff_R < 0) or (diff_L < 0):
            code = 0 # remove profile           
        else:
            code = 1 # keep profile
        
        # add info to gdf
        gdf.loc[(gdf['profielmet'] == name), 'filter_MV'] = code
    
   
        '''Apply a filter based on minimum profile height: it shouldn't be lower than -10 m NAP '''
        # Minimum profile height
        MV_min       = np.min(group['profiel'].values)
        
        # Apply filter based on height difference
        if MV_min < -10:
            code = 0 # remove profile           
        else:
            code = 1 # keep profile
        
        # add info to gdf
        gdf.loc[(gdf['profielmet'] == name), 'filter_MVmin'] = code


    # Apply filter 1
    gdf_sel  = gdf[gdf['filter_gaps'] == 1]
    profiles = np.unique(gdf['profielmet'].values)
    names    = np.unique(gdf[gdf['filter_gaps'] == 0]['profielmet'].values)
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because of missing AHN data at A24 or A28:\n')
        outfile.write(str(list(names)))
        outfile.write('\n')
    print(len(names), ' of ',len(profiles), ' profiles are dropped because of missing AHN data at A24 or A28. Check ', paths.txt_dropped, ' for the ID codes.')
    
    # Apply filter 2
    gdf      = gdf_sel.copy()
    gdf_sel  = gdf[gdf['filter_MV'] == 1]
    profiles = np.unique(gdf['profielmet'].values)
    names    = np.unique(gdf[gdf['filter_MV'] == 0]['profielmet'].values)
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because of AHN data at A24/A28 is lower than at A22L/A22R:\n')
        outfile.write(str(list(names)))
        outfile.write('\n')
    print(len(names), ' of ',len(profiles), ' profiles are dropped because of AHN data at A24/A28 is lower than at A22L/A22R. Check ', paths.txt_dropped, ' for the ID codes.')
        
    # Apply filter 3
    gdf      = gdf_sel.copy()
    gdf_sel  = gdf[gdf['filter_MVmin'] == 1]
    profiles = np.unique(gdf['profielmet'].values)
    names    = np.unique(gdf[gdf['filter_MVmin'] == 0]['profielmet'].values)
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because the minimum height is below -10 m NAP:\n')
        outfile.write(str(list(names)))
        outfile.write('\n')
    print(len(names), ' of ',len(profiles), ' profiles are dropped because the minimum height is below -10 m NAP. Check ', paths.txt_dropped, ' for the ID codes.')
    
    
    # remove remaining gaps
    gdf['profiel'] = np.where(gdf['profiel'] < -1e10, np.nan, gdf['profiel'])
    gdf['AHN']     = np.where(gdf['AHN'] < -1e10, np.nan, gdf['AHN'])
    
    # Save files
    gdf_sel.to_file(paths.shp_points_AHN_sel) 
     

def final_check(paths):
    # Input data
    gdf       = gpd.read_file(paths.shp_points_AHN_sel)
    df_HO     = gpd.read_file(paths.shp_HO)
    gdf['HO_ID'] = gdf['HO_ID'].astype(str)
        
    # Missing Hydro Objects
    HO_all  = np.unique(df_HO['CODE'].values)
    HO_prof = np.unique(gdf['HO_ID'].values)
    HO_code = [HO for HO in HO_all if HO not in HO_prof]
    print(len(HO_code), ' of ',len(HO_all), ' hydro objects have no profile.')
    
    # Missing Hydro Objects primary waterways
    df_HO_primair = df_HO[df_HO['CATEGORIEO']==1]
    HO_all  = np.unique(df_HO_primair['CODE'].values)
    HO_prof = np.unique(gdf['HO_ID'].values)
    HO_code = [HO for HO in HO_all if HO not in HO_prof]
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Primary hydro objects with no profile:\n')
        outfile.write(str(HO_code))
        outfile.write('\n')
    print(len(HO_code), ' of ',len(HO_all), ' primary hydro objects have no profile. Check ', paths.txt_dropped, ' for the ID codes.')
    
    
    df_HO_missing = df_HO.loc[df_HO['CODE'].isin(HO_code)]
    print(len(df_HO_missing[df_HO_missing ['LENGTE'] < 100]), ' primary hydro objects are shorter than 100 m.')
    