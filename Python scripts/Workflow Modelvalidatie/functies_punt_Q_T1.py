# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 13:14:34 2024

@author: PetraH
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import warnings
from shapely.geometry import Point, MultiLineString
import functies_punt_Q_T1_bokeh

warnings.filterwarnings("ignore")

t_min_coverage    = 3*365 # minimum number observation points required to validate
vmin              = 0
vmax              = 200

def getLineCoords(row, geom, coord_type):
    if isinstance(row[geom], MultiLineString):
        empty_l = []
        return empty_l
    else:
        if coord_type == 'x':
            return list( row[geom].coords.xy[0] )
        elif coord_type == 'y':
            return list( row[geom].coords.xy[1] )
    
def fun_remove(df, name):
    df_sel = df.copy()
    for xy, group in df.groupby(['X','Y']):
        if (name in group['TYPE'].values) & (len(group['TYPE'].values)>0):
            ind = np.where((df_sel['X']==xy[0]) & (df_sel['Y']==xy[1]) & (df_sel['TYPE']!=name))[0]
            ind = df_sel.iloc[ind].index.values            
            df_sel = df_sel.drop(ind)
            
    return df_sel

def add_stats(df_sites_Q, df_data_Q):
    
    # selecteer data
    for i in range(0,len(df_sites_Q)):            
        # selecteer locaties die inbegrepen zijn in de data csv
        ID = df_sites_Q.iloc[i]['ID']
        
        # extract data & add to xds
        df_sel = df_data_Q[np.append(['YEAR'], ID)]
        df_sel = df_sel.astype({ID: 'float'})                               
        df_sel[ID] = df_sel[ID].replace(-999,np.nan)
        df_sel = df_sel.dropna()
        
        # Estimate stats
        if len(df_sel) > t_min_coverage:
            df_selmax   = df_sel.groupby('YEAR', dropna=True).max()
            data_T1     = float(df_selmax.mean().values)
        else: 
            data_T1 = -999
            
        # add to dataframe
        df_sites_Q.loc[(df_sites_Q['ID'] == ID), 'Data_T1'] = data_T1
    
    return df_sites_Q

def Sobek_2_shp_Q(df_Q, shp_reach):
    # Link Sobek data to shapefile: Debiet
    rmax        = len(shp_reach['ID'].values)
    shp_reach_Q = np.zeros((rmax,))*np.nan
    for i in range(0,len(shp_reach_Q)):
        ind = np.where((shp_reach['ID'].values[i]==df_Q['ID']))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_reach_Q[i] = float(df_Q['Model_T1'].values[ind])
    shp_reach['Model_T1'] = np.round(shp_reach_Q,1)
    shp_reach = shp_reach[shp_reach['Model_T1'].notna()]
    
    return shp_reach

def fun_remove_duplicates(gdf_data):
    # Step 1 - collect all rows that are *not* duplicates (based on ID)
    non_duplicates_to_keep = gdf_data.drop_duplicates(subset='ID_left', keep=False)

    # Step 2a - identify *all* rows that have duplicates (based on ID, keep all)
    sub_df = gdf_data[gdf_data.duplicated('ID_left', keep=False)]

    # Step 2b - of those duplicates, discard all that have "0" in column
    duplicates_to_keep = sub_df[sub_df['Model_T1'] != 0]

    # join the 2 sets
    gdf_data = pd.concat([non_duplicates_to_keep, duplicates_to_keep])
    
    # sort
    gdf_data = gdf_data.sort_values(by=['ID_left'], ascending=True)
    
    return gdf_data
    
def main(paths, months, years):
    ''''FEWS-WIS data'''
    # Input data
    df_sites_Q = pd.read_csv(paths.df_Q_WISsites)
    df_sites_Q = df_sites_Q.rename(columns={"LOC_ID": "ID"}).drop(0)
    df_data_Q  = pd.read_csv(paths.df_Q_WISdata, delimiter=',',low_memory=False)    
    df_data_Q  = df_data_Q.rename(columns={"GMT+1": "time"}).drop(0)
    df_data_Q  = df_data_Q.astype({'time':'datetime64[ns]'}) 
    df_data_Q  = df_data_Q.set_index('time').astype('float')
    df_data_Q[df_data_Q<-900] = -999
    df_data_Q  = df_data_Q.replace(-999,np.nan).reset_index()
    df_data_Q['MONTH'] = np.array([pd.to_datetime(t).month for t in df_data_Q.time.values])
    df_data_Q['YEAR'] = np.array([pd.to_datetime(t).year for t in df_data_Q.time.values])    
    shp_insteek= gpd.read_file(paths.shp_insteek)
    shp_afvoer = gpd.read_file(paths.shp_afvoer)
    
    # Select type: debietmeter, pompvijzel
    df_sites_Q = df_sites_Q.loc[(df_sites_Q['TYPE']=='debietmeter') | (df_sites_Q['TYPE']=='pompvijzel')]
    
    # Select sites with data
    idx_data   = [i for i in range(0,len(df_sites_Q)) if df_sites_Q.ID.values[i] in df_data_Q.columns]
    df_sites_Q = df_sites_Q.iloc[idx_data] 
    
    # bereken statistiek
    df_data_Q_sel  = df_data_Q.copy()
    df_data_Q_sel  = df_data_Q_sel[df_data_Q_sel['MONTH'].isin(months)]
    df_data_Q_sel  = df_data_Q_sel[df_data_Q_sel['YEAR'].isin(years)]
    data_T1  = add_stats(df_sites_Q, df_data_Q_sel)
    
    # Create geodataframe
    cols         = ['ID','LOC_NAME','X','Y','TYPE','Data_T1','geometry']    
    gdf_data     = gpd.GeoDataFrame(data_T1, geometry=[Point(xy) for xy in zip(data_T1.X,data_T1.Y)])
    gdf_data     = gdf_data[cols]
    gdf_data.crs ='EPSG:28992'
    
    
    # Clip op primary waterways
    shp_insteek = shp_insteek.loc[shp_insteek['CATEGORIEO']==1]
    gdf_data    = gpd.clip(gdf_data, shp_insteek)
    
    # Save data
    gdf_data.to_file(paths.gdf_data)
    
    '''Model data'''
    # data inlezen
    shp_reach    = gpd.read_file(paths.shp_reach)
    df_Q_model   = pd.read_csv(paths.df_Q_model,skiprows=4, sep = ',')
    df_Q_model.columns = ['ID','Model_T1']
    df_Q_model   = df_Q_model.astype({'ID':'str','Model_T1':'float'}) 
    
    # Spaties in kolom namen verwijderen
    shp_reach.columns = [col.strip() for col in shp_reach.columns]
    
    # Link Sobek data aan shapefile
    gdf_model = Sobek_2_shp_Q(df_Q_model, shp_reach)
    
    # Save file
    gdf_model.crs='EPSG:28992'
    gdf_model.to_file(paths.gdf_model)
    
    '''Data vs. Model'''
    gdf_data = gpd.read_file(paths.gdf_data)
    gdf_model = gpd.read_file(paths.gdf_model)
    
    gdf_joined = gpd.sjoin_nearest(gdf_data, gdf_model, distance_col="distances")
    gdf_joined['dQ'] = gdf_joined['Model_T1']/gdf_joined['Data_T1']*100
    gdf_joined['dQ'] = np.where((gdf_joined['Model_T1'] == -999) | (gdf_joined['Data_T1'] == -999), -999, gdf_joined['dQ'])  
    gdf_joined['dQ'] = np.where(gdf_joined['Data_T1']==0, -999, gdf_joined['dQ'])            
    gdf_joined['dQ'] = np.round(gdf_joined['dQ'])
    gdf_joined       = gdf_joined[gdf_joined['dQ']!=-999]
    
    # Remove duplicates (remove row with Q=0)
    gdf_joined = fun_remove_duplicates(gdf_joined)
    
    # Save file
    gdf_joined.to_file(paths.gdf_joined)
    
    # spatial join (nodes & shapefile)
    shp_node_joined = gpd.sjoin(gdf_joined[['dQ','geometry']], shp_afvoer[['CODE','geometry']], how="left")
    shp_node_joined['dQabs'] = abs(shp_node_joined['dQ'])
    shp_node_joined = shp_node_joined.loc[shp_node_joined.groupby('CODE')['dQabs'].idxmax().values]
    shp_afvoer      = shp_afvoer.merge(shp_node_joined[['CODE','dQ']], on='CODE', how='left')
    
    # Data plotten
    # bokeh_fig(paths, shp_afvoer, gdf_joined)
    functies_punt_Q_T1_bokeh.bokeh_fig(paths, shp_afvoer, gdf_joined, df_data_Q, vmin, vmax)
    
    
    