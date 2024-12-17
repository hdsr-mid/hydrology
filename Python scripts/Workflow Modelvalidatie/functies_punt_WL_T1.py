# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 13:14:34 2024

@author: PetraH
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
import warnings
from shapely.geometry import Point, MultiLineString
warnings.filterwarnings("ignore")

import functies_punt_WL_T1_bokeh

max_dist_2_struc  = 50 # max distance to pump/weir/culvert -> used for conditional approach in code
t_min_coverage    = 3*365 # minimum number observation points required to validate
vmin              = -1
vmax              = 1

    
def fun_stat_data(df_sites, df_data):
    
    # selecteer data
    for i in range(0,len(df_sites)):            
        # selecteer locaties die inbegrepen zijn in de data csv
        ID = df_sites.iloc[i]['ID']
        
        # extract data & add to xds
        df_sel = df_data[np.append(['YEAR'], ID)]
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
        df_sites.loc[(df_sites['ID'] == ID), 'Data_T1'] = data_T1
    
    return df_sites

        
def getLineCoords(row, geom, coord_type):
    if isinstance(row[geom], MultiLineString):
        empty_l = []
        return empty_l
    else:
        if coord_type == 'x':
            return list( row[geom].coords.xy[0] )
        elif coord_type == 'y':
            return list( row[geom].coords.xy[1] )
   
def custom_resampler(arraylike):
    return np.nanmean(arraylike)
    
def main(paths, months, years):
    
    ''''FEWS-WIS data'''
    # Input data
    df_sites = pd.read_csv(paths.df_WL_WISsites)   
    df_sites = df_sites.rename(columns={"LOC_ID": "ID"}).drop(0)
    df_data  = pd.read_csv(paths.df_WL_WISdata, delimiter=',',low_memory=False)    
    df_data  = df_data.rename(columns={"GMT+1": "time"}).drop(0)
    cols     = [c for c in df_data.columns if c != 'time']
    df_data  = df_data.astype({'time':'datetime64[ns]'})   
    df_data[cols]    = df_data[cols].astype(np.float64)
    df_data['MONTH'] = np.array([pd.to_datetime(t).month for t in df_data.time.values])
    df_data['YEAR']  = np.array([pd.to_datetime(t).year for t in df_data.time.values])
    
    # df_data    = df_data.set_index('time').resample('D').apply(custom_resampler).reset_index()
    shp_insteek= gpd.read_file(paths.shp_insteek).to_crs(28992)
    shp_afvoer = gpd.read_file(paths.shp_afvoer)
    
    # Select sites with data
    idx_data   = [i for i in range(0,len(df_sites)) if df_sites.ID.values[i] in df_data.columns]
    df_sites = df_sites.iloc[idx_data] 
    
    # bereken statistiek
    df_data_sel  = df_data.copy()
    df_data_sel  = df_data_sel[df_data_sel['MONTH'].isin(months)]
    df_data_sel  = df_data_sel[df_data_sel['YEAR'].isin(years)]
    data_stat    = fun_stat_data(df_sites, df_data_sel)
    
    # omzetten naar shapefile
    gdf_data  = gpd.GeoDataFrame(data_stat, geometry=[Point(xy) for xy in zip(data_stat.X,data_stat.Y)])
    gdf_data.crs='EPSG:28992'
    
    # clip op primaire watergangen    
    shp_insteek = shp_insteek.loc[shp_insteek['CATEGORIEO']==1]
    gdf_data = gpd.clip(gdf_data, shp_insteek)
    
    # insteek code toevoegen via spatial join
    gdf_data = gpd.sjoin(gdf_data, shp_insteek[['CODE','geometry']].rename(columns={"CODE": "IS_CODE"}), how="left")
    
    # opslaan als shapefile
    gdf_data.to_file(paths.gdf_data)
    
    '''Model data'''
    # data inlezen
    cols               = ['ID','TYPE','X','Y','geometry']
    shp_node           = gpd.read_file(paths.shp_node)[cols]
    shp_node_peilgrens = shp_node[(shp_node['TYPE']=='Flow - Weir') | (shp_node['TYPE']=='Flow - Culvert') | (shp_node['TYPE']=='Flow - Pump Station')]
    shp_node           = gpd.clip(shp_node, shp_insteek) # met name de storage nodes worden hiermee verwijderd
    df_WL              = pd.read_csv(paths.df_model,skiprows=4, delimiter = ',').reset_index()
    df_WL.columns      = ['index','ID','WL']    
    
    # Link Sobek data to shapefile: Waterniveau
    shp_node=shp_node[~shp_node["ID"].str.contains('stor', na=False)] # remove stor-nodes
    nmax         = len(shp_node['ID'].values)
    shp_node_WL  = np.zeros((nmax, 1))-999
    for i in range(0,len(df_WL)):
        ind = np.where((df_WL['ID'][i]==shp_node['ID'].values))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_node_WL[ind]  = df_WL['WL'].values[i]         
    shp_node['Model_T1']  = np.round(shp_node_WL,1)
    shp_node              = shp_node[shp_node['Model_T1']!=-999]
    
    # opslaan als shapefile
    shp_node.to_file(paths.gdf_model)    
    
    '''Data vs. Model'''
    gdf_data  = gpd.read_file(paths.gdf_data).to_crs(28992)
    gdf_model = gpd.read_file(paths.gdf_model).to_crs(28992)
    
    # join: model - data
    gdf_joined = gpd.GeoDataFrame()
    
    for gg in tqdm(range(0,len(gdf_data))):
        row = gdf_data[gg:gg+1]        
        distance = gpd.sjoin_nearest(row, shp_node_peilgrens, distance_col="distances")['distances'].min()        
        if distance < max_dist_2_struc:
            # vind bijbehorend insteekvlak
            shp_insteek_gg = shp_insteek.loc[shp_insteek['CODE'].values==row['IS_CODE'].values][['CODE','geometry']]
            # voeg kleine buffer toe om geometrische afrondingen te omzeilen
            shp_insteek_gg['geometry'] = shp_insteek_gg.geometry.buffer(0.01)
            # vind bijbehorende modelpunten in insteekvlak
            gdf_model_gg = gpd.clip(gdf_model, shp_insteek_gg)            
            # vind dichtstbijzijnde modelpunt uit bovenstaandje collectie
            joined_gg = gpd.sjoin_nearest(row, gdf_model_gg, distance_col="distances")            
        else:
            # vind dichtstbijzijnde modelpunt uit gehele collectie
            joined_gg = gpd.sjoin_nearest(row, gdf_model, distance_col="distances")
        gdf_joined= pd.concat([gdf_joined,joined_gg], ignore_index=True)
    gdf_joined['dWL'] = gdf_joined['Model_T1'] - gdf_joined['Data_T1']   
    gdf_joined['dWL'] = np.where((gdf_joined['Model_T1'] == -999) | (gdf_joined['Data_T1'] == -999), -999, gdf_joined['dWL'])        
    gdf_joined        = gdf_joined[gdf_joined['dWL']!=-999]
    gdf_joined.to_file(paths.gdf_joined)
    
    # spatial join (nodes & shapefile)
    shp_node_joined = gpd.sjoin(gdf_joined[['dWL','geometry']], shp_afvoer[['CODE','geometry']], how="left")
    shp_node_joined['dWLabs'] = abs(shp_node_joined['dWL'])
    shp_node_joined = shp_node_joined.loc[shp_node_joined.groupby('CODE')['dWLabs'].idxmax().values]
    shp_afvoer      = shp_afvoer.merge(shp_node_joined[['CODE','dWL']], on='CODE', how='left')
    
    # Data plotten
    # shp_regios = gpd.read_file(paths.shp_regios)
    # region1    = shp_regios[shp_regios['REGIO_WSB'] == 'KROMME RIJN']
    # region2    = shp_regios[shp_regios['REGIO_WSB'] != 'KROMME RIJN']

    # functies_punt_WL_T1_bokeh.bokeh_fig(paths, shp_afvoer, gdf_joined, region1, df_data, vmin, vmax, 'KROMME RIJN')
    # functies_punt_WL_T1_bokeh.bokeh_fig(paths, shp_afvoer, gdf_joined, region2, df_data, vmin, vmax, 'OVERIGE')
    functies_punt_WL_T1_bokeh.bokeh_fig(paths, shp_afvoer, gdf_joined, df_data, vmin, vmax, 'all')
    