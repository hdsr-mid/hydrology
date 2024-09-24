# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:47:43 2024

@author: PetraH
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

class paths():
    # Input data
    root        = os.getcwd()
    shp_reach   = os.path.join(root,"Input","RchSegments.shp")
    shp_node    = os.path.join(root,"Input",'3b_nod.shp')
    shp_peil    = os.path.join(root,'Input','BR-VS_Peilgebieden.shp')
    df_Q        = os.path.join(root,'Input','Q_RUPROF_150724.csv')
    df_WL       = os.path.join(root,'Input','WL_RUPROF_150724.csv')
    
    # Output data
    shp_out    = os.path.join(root,'Output','resultaten.gpkg')
    

def Sobek_2_shp_Q(df_Q, shp_reach):
    # Link Sobek data to shapefile: Debiet
    df_Q.columns= ['ID','Q']
    rmax        = len(shp_reach['ID'].values)
    shp_reach_Q = np.zeros((rmax,))*np.nan
    lwidth      = np.zeros((rmax,))*np.nan
    Q_upper     = np.percentile(df_Q.Q.values,95)
    Q_lower     = np.percentile(df_Q.Q.values,5)
    upper_limit = 5
    for i in tqdm(range(0,len(shp_reach_Q))):
        ind = np.where((shp_reach['ID'].values[i]==df_Q['ID']))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_reach_Q[i] = float(df_Q['Q'].values[ind])
            l_width        = (float(shp_reach_Q[i]) - Q_lower)/(Q_upper - Q_lower)
            lwidth[i]      = np.nanmin([np.nanmax([0, l_width  * upper_limit]),upper_limit]) + 0.5
    shp_reach['Q'] = shp_reach_Q
    shp_reach['lw']= lwidth
    
    return shp_reach, Q_upper, Q_lower


def Sobek_2_shp_WL(df_WL, shp_node):
    # Link Sobek data to shapefile: Waterlevel
    df_WL.columns  = ['ID','WL']
    
    nmax           = len(shp_node['ID'].values)
    shp_node_WL    = np.zeros((nmax, 1))*np.nan
    for i in tqdm(range(0,len(shp_node_WL))):
        ind = np.where((shp_node['ID'].values[i]==df_WL['ID']))[0]
        if len(ind)>0: 
            ind               = int(ind[0]) #Attention: soms zijn er meerdere waarden voor 1 knoop; voor nu is de 1e waarde genomen
            shp_node_WL[i]    = df_WL['WL'].values[ind]                
    shp_node['WL']     = shp_node_WL
       
    return shp_node


def point_to_poly(shp_node, shp_poly):
    CODE = 'ID'
    
    # Add max(WLmax) to shp_poly
    WL_max_peil  = np.zeros(len(shp_poly['CODE']),)-99
    for i in tqdm(range(0,len(shp_poly['CODE']))):
        idx = np.where(shp_node.within(shp_poly.iloc[i].geometry))[0]
        
        if len(idx)>0:
            WL_max_peil[i] = np.nanmax(shp_node.iloc[idx]['WL'].values)
    shp_poly['WLmax']  = WL_max_peil
    
    return shp_poly

if __name__ == "__main__":
    
    # Input data
    shp_peil     = gpd.read_file(paths.shp_peil)
    shp_reach    = gpd.read_file(paths.shp_reach)
    shp_node     = gpd.read_file(paths.shp_node)
    df_Q         = pd.read_csv(paths.df_Q,skiprows=4)
    df_WL        = pd.read_csv(paths.df_WL,skiprows=4)
    
    # Edit shapefile column notation
    columns_shp_reach = [c.replace('  ','').replace('  ','').replace(' ','') for c in shp_reach.columns]
    columns_shp_node  = [c.replace('  ','').replace('  ','').replace(' ','') for c in shp_node.columns]
    shp_reach.columns = columns_shp_reach
    shp_node.columns  = columns_shp_node
    
    # Link Sobek data to shapefile: Debiet
    [shp_reach,Q_upper, Q_lower] = Sobek_2_shp_Q(df_Q, shp_reach)
    
    # Link Sobek data to shapefile: Waterlevel
    shp_node = Sobek_2_shp_WL(df_WL, shp_node)
    
        
    # Add max(WL) value to shp_poly
    shp_peil = point_to_poly(shp_node, shp_peil)
    
    # save files
    shp_node.to_file(paths.shp_out, driver="GPKG",layer = 'Sobek_WL')
    shp_reach.to_file(paths.shp_out, driver="GPKG",layer = 'Sobek_Q')
    shp_peil.to_file(paths.shp_out, driver="GPKG",layer = 'Peilgebieden')
    
    
    