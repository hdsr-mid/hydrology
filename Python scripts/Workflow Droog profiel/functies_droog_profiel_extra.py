# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:26:20 2024

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

def format_aanpassingen(paths):
    '''Eit tile input file if needed
    In this case, the codes A22L en A22R were added (before it was called A22)
    Also, the column "afstand" was added (distance)
    
    Maak kleine aanpassingen aan de format
    - schijdingsteken ";"
    - decimaal teken ","
    - kolomnamen: profielmeting, punttype, afstand, slibhoogte, X, Y, vastebodem
    - kolom "punttype": A23, A24, A22L, A99, A22R, A28, A29
    
    In deze script zijn de volgende aanpassingen gemaakt aan de input data:
        - punttype: van integer nummer naar A23 - A29
        - afstand : berekend op basis van XY coordinaten, met als start (afstand = 0): helemaal links    
    
    '''
    
    df = pd.read_csv(paths.txt_profielen_temp,delimiter=',')

    # Edit 1: column "punttype"
    df['punttype'] = df['punttype'].astype(str)
    
    for profielmet, group in df.groupby('profielmeting'): 
        punttype = ['A' + i for i in group.punttype.values]
        ind = np.where(np.array(punttype)=='A22')[0]
        punttype[ind[0]] = 'A22L'
        punttype[ind[1]] = 'A22R'
        df.loc[df['profielmeting']==profielmet,'punttype'] = punttype
    
    # Add geometry (needed for next step)
    gdf  = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.X,df.Y)])
    gdf.crs='EPSG:28992'
    
    # Edit 2: extra column "afstand"
    df['afstand'] = 0
    for profielmet, group in gdf.groupby('profielmeting'): 
        reference_point = group[group['punttype'] == group['punttype'].values[0]].geometry.iloc[0]
        
        afstand = group.geometry.apply(lambda x: x.distance(reference_point))
        afstand = np.round(afstand,2)
        df.loc[df['profielmeting']==profielmet,'afstand'] = afstand
    
    df.to_csv(paths.txt_profielen)
    
def selectie_nat_profiel(paths):
    ''' Select the wet part of the profile (to mimick WIT data)'''
    
    # Input data
    df   = pd.read_csv(paths.txt_profielen)
    gdf  = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.X,df.Y)])
    gdf.crs='EPSG:28992'
    gdf.to_file(paths.shp_profielen_gemeten)
    
    # Select wet points
    df_nat_profiel    = df[(df['punttype']=='A99') | (df['punttype']=='A22L') | (df['punttype']=='A22R')]
    
    # Save as shapefile
    gdf  = gpd.GeoDataFrame(df_nat_profiel, geometry=[Point(xy) for xy in zip(df_nat_profiel.X,df_nat_profiel.Y)])
    gdf.crs='EPSG:28992'
    gdf.to_file(paths.shp_profielen_nat)
def extra_afstd_droog(paths):
    gdf = gpd.read_file(paths.shp_points_incl_AHN)
    gdf.crs='EPSG:28992'
    shp = gpd.read_file(paths.shp_water)
    shp.crs='EPSG:28992'
    
    # Groepeer de punten op basis van de kolom 'profielmet'
    grouped = gdf.groupby('profielmet')
    gdf_L = gdf[gdf['punttype'] == 'A22L']
    gdf_R = gdf[gdf['punttype'] == 'A22R']
    
    # Loop door elke profielmet-code
    for name, group in grouped:
        
        # links
        # afstand insteek vs. waterlijn van gemeten punten
        p_insteek    = group[group['punttype'] == 'A24'].geometry.iloc[0]
        p_water      = group[group['punttype'] == 'A22L'].geometry.iloc[0]
        afstd        = p_insteek.distance(p_water)
        gdf_L.loc[(gdf_L['profielmet'] == name), 'A24_A22L'] = afstd
        # afstand naar waterlijn (gemeten punt A22 vs. BR waterlijn)
        afstd        = shp.boundary.distance(p_water).sort_values().values[0]
        gdf_L.loc[(gdf_L['profielmet'] == name), 'A22L_BR'] = afstd
        
        
        # rechts
        # afstand insteek vs. waterlijn van gemeten punten
        p_insteek    = group[group['punttype'] == 'A28'].geometry.iloc[0]
        p_water      = group[group['punttype'] == 'A22R'].geometry.iloc[0]
        afstd        = p_insteek.distance(p_water)
        gdf_R.loc[(gdf_R['profielmet'] == name), 'A28_A22R'] = afstd
        # afstand naar waterlijn (gemeten punt A22 vs. BR waterlijn)
        afstd        = shp.boundary.distance(p_water).sort_values().values[0]
        gdf_R.loc[(gdf_R['profielmet'] == name), 'A22R_BR'] = afstd
        
        
    gdf_prof = gpd.GeoDataFrame(pd.concat([gdf_L, gdf_R], ignore_index=True))    
    gdf_prof.to_file(paths.shp_points_dist)    

def add_HO_ID(paths, gdf):
    ''' Add the code and smallest distance to a Hydro object to the profiles'''
    
    # Input data
    shp_HO = gpd.read_file(paths.shp_HydroObject)
    shp_HO = shp_HO[['CODE','CATEGORIEO','geometry']]
    
    # Group points based on the column 'profielmet'
    grouped = gdf.groupby('profielmet')
    
    # Loop through each profielmet-code
    for name, group in grouped:    
        # Spatial join
        gdf_joined = gpd.sjoin_nearest(group, shp_HO, distance_col="distances")
        
        # get smallest distance
        distance = np.min(gdf_joined.distances.values)
        
        # get code of Hydro object
        CODE = gdf_joined.CODE.values[0]
        
        # add info to gdf
        gdf.loc[(gdf['profielmet'] == name), 'HO_distance'] = distance
        gdf.loc[(gdf['profielmet'] == name), 'HO_CODE'] = CODE
    
    return gdf  
    
def fun_filter_HO(paths):
    '''Apply a filter based on the variability (of the height at the profile edge; insteek) within a Hydro object'''    
    
    # Input data
    gdf    = gpd.read_file(paths.shp_points_AHN_sel)
    
    # Add Hydro object code to profile
    add_HO_ID(paths, gdf)
    
    # Group points based on the column 'profielmet'
    grouped = gdf.groupby('HO_CODE')
    
    # Get list of profiles (one row per profile)
    profiles = gdf[gdf['punttype'] == 'A24']
    
    # Loop through Hydro objects
    for name, group in grouped:  
        # Criteria: difference in height at profile edge (insteek hoogte) left side
        p_insteek    = group[group['punttype'] == 'A24']
        diff_L       = np.max(p_insteek['profiel'].values) - np.min(p_insteek['profiel'].values)
        
        # Criteria: difference in height at profile edge (insteek hoogte) right side
        p_insteek    = group[group['punttype'] == 'A28']
        diff_R       = np.max(p_insteek['profiel'].values) - np.min(p_insteek['profiel'].values)
        
        # Criteria: difference in maximum profile height
        grouped_profiles = group.groupby('profielmet')
        for profile_name, profile in grouped_profiles:
            height  = np.max(profile['profiel'].values) - np.min(profile['profiel'].values)
            group.loc[(group['profielmet'] == profile_name), 'height'] = height
        diff_height = np.max(group['height'].values) - np.min(group['height'].values)
        
        # Apply filter based on height difference
        if (diff_R > 0.5) or (diff_L > 0.5) or (diff_height > 0.5):
            print('HO')
            print(name)
            print(group[['profielmet', 'HO_CODE','punttype', 'profiel','height']])
            print('')
            code = 0 # remove profile
        else:
            code = 1 # keep profile
        
        
        # add info to gdf
        for profile_ID in group['profielmet'].values:
            gdf.loc[(gdf['profielmet'] == profile_ID), 'filter_HO'] = code
        
    # Check
    profiles = gdf[gdf['punttype'] == 'A24']
    print('HO',np.sum(profiles.filter_HO.values), len(profiles))
    
    # Save file
    if os.path.exists(paths.shp_points_AHN_sel):
        os.remove(paths.shp_points_AHN_sel)
        gdf.to_file(paths.shp_points_AHN_sel)        
    else:
        gdf.to_file(paths.shp_points_AHN_sel)      