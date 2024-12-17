# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 13:14:34 2024

@author: PetraH
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import warnings
from shapely.geometry import Point
from pathlib import Path
import hkvsobekpy as hkv
import functies_punt_Qseries_bokeh
warnings.filterwarnings("ignore")

t_min_coverage    = 3*365 # minimum number observation points required to validate

def add_stats(gdf_joined, df_sites_Q, df_data_Q, sub_data):
    # Create empty columns
    gdf_joined['NS']     = -999
    gdf_joined['RMSE']   = -999
    gdf_joined['Q_MOD']  = -999
    gdf_joined['Q_OBS']  = -999
    gdf_joined['dQ']     = -999
    
    for i in range(0,len(gdf_joined)):    
        ID_data  = gdf_joined.iloc[i]['ID_left']
        ID_model = gdf_joined.iloc[i]['ID_right']
        df_data = df_data_Q[['time',ID_data]]
        df_model= sub_data[ID_model].reset_index()
        df_model.columns =['time', ID_model]
        df_model = df_model.set_index('time').resample('D').mean().reset_index()
        df_data  = df_data.set_index('time').resample('D').mean().reset_index()
        df_sel   = df_model.merge(df_data, on='time', how='left')
        
        # Performance
        OBS    = df_sel.dropna()[ID_data].values
        MOD    = df_sel.dropna()[ID_model].values
        OBS    = OBS[np.isnan(OBS)==0]
        MOD    = MOD[np.isnan(OBS)==0]
        OBSavg = np.nanmean(df_sel[ID_data].dropna().values)
        MODavg = np.nanmean(df_sel[ID_model].dropna().values)
        if np.isnan(OBSavg): OBSavg = -999
        
        gdf_joined.loc[(gdf_joined['ID_left'] == ID_data), 'Q_MOD']  = MODavg
        gdf_joined.loc[(gdf_joined['ID_left'] == ID_data), 'Q_OBS']  = OBSavg
        
            
        if len(OBS) > t_min_coverage:        
            RMSE = np.round(np.sqrt(np.sum((OBS - MOD)**2)/len(OBS)),2)
            NS   = np.round(1 - np.sum((OBS-MOD)**2)/np.sum((OBS-np.mean(OBS))**2),2)
            if np.isinf(NS): NS = np.nan
            
            dQ     = MODavg/OBSavg*100
            if (OBSavg == -999) | (OBSavg == 0): dQ = -999
            
            gdf_joined.loc[(gdf_joined['ID_left'] == ID_data), 'NS']     = NS
            gdf_joined.loc[(gdf_joined['ID_left'] == ID_data), 'RMSE']   = RMSE
            gdf_joined.loc[(gdf_joined['ID_left'] == ID_data), 'dQ']     = dQ
            
    return gdf_joined

def custom_resampler(arraylike):
    return np.nanmean(arraylike)
 
# spatial join (nodes & shapefile)
def node_2_poly(metric,gdf_joined,shp_afvoer):
     gdf_joined_sel           = gdf_joined[gdf_joined[metric]!=-999]
     shp_node_joined          = gpd.sjoin(gdf_joined_sel[[metric,'geometry']], shp_afvoer[['CODE','geometry']], how="left")
     shp_node_joined['abs']   = abs(shp_node_joined[metric])
     shp_node_joined          = shp_node_joined.dropna()
     shp_node_joined          = shp_node_joined.loc[shp_node_joined.groupby('CODE')['abs'].idxmax().values]
     shp_afvoer               = shp_afvoer.merge(shp_node_joined[['CODE',metric]], on='CODE', how='left')
     return shp_afvoer   


def main(paths, months, years):
    ''''FEWS-WIS data'''
    # Input data
    df_sites_Q = pd.read_csv(paths.df_Q_WISsites)    
    df_sites_Q = df_sites_Q.rename(columns={"LOC_ID": "ID"}).drop(0)
    df_data_Q  = pd.read_csv(paths.df_Q_WISdata, delimiter=',',low_memory=False)
    df_data_Q  = df_data_Q.rename(columns={"GMT+1": "time"})
    df_data_Q  = df_data_Q.astype({'time':'datetime64[ns]'})  
    df_data_Q  = df_data_Q.set_index('time').astype('float')
    df_data_Q[df_data_Q<-900] = -999
    df_data_Q  = df_data_Q.replace(-999,np.nan).reset_index()    
    df_data_Q['MONTH'] = np.array([pd.to_datetime(t).month for t in df_data_Q.time.values])
    df_data_Q['YEAR'] = np.array([pd.to_datetime(t).year for t in df_data_Q.time.values])
    df_data_Q  = df_data_Q[df_data_Q['MONTH'].isin(months)]
    df_data_Q  = df_data_Q[df_data_Q['YEAR'].isin(years)]
    df_data_Q  = df_data_Q.set_index('time').resample('D').apply(custom_resampler).reset_index()
    shp_insteek= gpd.read_file(paths.shp_insteek)
    shp_afvoer = gpd.read_file(paths.shp_afvoer)
    
    # Select type: debietmeter, pompvijzel
    df_sites_Q = df_sites_Q.loc[(df_sites_Q['TYPE']=='debietmeter') | (df_sites_Q['TYPE']=='pompvijzel')]
    idx_data   = [i for i in range(0,len(df_sites_Q)) if df_sites_Q.ID.values[i] in df_data_Q.columns]
    df_sites_Q = df_sites_Q.iloc[idx_data]    
    
    # Create geodataframe
    gdf_data  = gpd.GeoDataFrame(df_sites_Q, geometry=[Point(xy) for xy in zip(df_sites_Q.X,df_sites_Q.Y)])
    gdf_data.crs='EPSG:28992'
    
    # Clip op primaire watergangen
    shp_insteek = shp_insteek.loc[shp_insteek['CATEGORIEO']==1].buffer(1)
    gdf_data = gpd.clip(gdf_data, shp_insteek)
    
    # Save file
    cols     = ['ID','LOC_NAME','X','Y','TYPE','IRIS_ID','geometry']
    gdf_data = gdf_data[cols]
    gdf_data.to_file(paths.gdf_data)
    
    '''Model data'''
    # data inlezen
    shp_reach    = gpd.read_file(paths.shp_reach)
    shp_reach.columns = [col.strip() for col in shp_reach.columns]
    # shp_reach.index = shp_reach['ID']
    shp_reach.crs='EPSG:28992'
    shp_reach.drop(['NAME', "TYPE", 'PARENTID', 'ID_FROM', 'ID_TO', 'USERID'], axis=1, inplace=True)
    
    cases        = [1,2,3,4,5,6,7,8,9,10,11,12]#
    subhislist = []
    for case in cases:
        pad_to_his = str(Path(paths.data_sobek) / str(case) / 'reachseg.his')
        reachseg = hkv.read_his.ReadMetadata(pad_to_his)    
        reach_segments = reachseg.DataFrame()
        reach_segments = reach_segments.iloc[32:,:]
        locs= reachseg.GetLocations()
        params= reachseg.GetParameters()
        reach_segments = reach_segments.loc[:,params[params.index('Discharge mean(m³/s)')]]
        reach_segments.head()
        subhislist.append(reach_segments)
    all_data = pd.concat(subhislist)
    all_data['MONTH'] = all_data.index.month
    all_data['YEAR']  = all_data.index.year
    sub_data  = all_data[all_data['MONTH'].isin(months)]
    sub_data  = all_data[all_data['YEAR'].isin(years)]
    idx_model = [i for i in range(0,len(shp_reach)) if shp_reach.ID.values[i] in sub_data.columns]
    gdf_model = shp_reach.iloc[idx_model]
    
    # Save files
    gdf_model = gdf_model.reset_index()
    gdf_model.to_file(paths.gdf_model)
    
    '''Data vs. Model'''
    # Estimate performance metric
    gdf_joined = gpd.sjoin_nearest(gdf_data, gdf_model, distance_col="distances")
    gdf_joined = add_stats(gdf_joined, df_sites_Q, df_data_Q, sub_data)
    
    # Save file
    gdf_joined.to_file(paths.gdf_joined)
    
    gdf_data = gpd.read_file(paths.gdf_data)
    gdf_model = gpd.read_file(paths.gdf_model)
    gdf_joined = gpd.read_file(paths.gdf_joined)
    
    # Kaart plotten
    metric    = 'dQ'
    vmin       = 0
    vmax       = 200
    
    shp_afvoer_joined    = node_2_poly(metric,gdf_joined,shp_afvoer)
         
    functies_punt_Qseries_bokeh.main(paths, shp_afvoer_joined, gdf_joined, sub_data, df_data_Q, metric, vmin, vmax)
        
        