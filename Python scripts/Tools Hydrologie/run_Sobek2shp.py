# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:47:43 2024

@author: PetraH
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import warnings
warnings.filterwarnings("ignore")

class paths():
    # Input data
    root        = os.getcwd()
    shp_reach   = os.path.join(root,"Input","RchSegments.shp")
    shp_node    = os.path.join(root,"Input",'3b_nod.shp')
    shp_peil    = os.path.join(root,'Input','BR-VS_Peilgebieden_clip.shp')
    df_Q        = os.path.join(root,'Input','Q.txt')
    df_WL       = os.path.join(root,'Input','WL.txt')
    
    # Output data
    Sobek_WL   = os.path.join(root,'Output','Sobek_WL.shp')
    Sobek_Q    = os.path.join(root,'Output',"Sobek_Q.shp")
    shp_peil_WL= os.path.join(root,'Output','BR-VS_Peilgebieden_clip_WL.shp')
    

def Sobek_2_shp_Q(df_Q, shp_reach):
    # Link Sobek data to shapefile: Debiet
    df_Q        = df_Q.drop(['Unnamed: 0','Location:'], axis=1)
    Q           = df_Q.iloc[0].reset_index() # timestep 0 -> change accordingly
    Q.columns   = ['ID','Q']
    rmax        = len(shp_reach['ID'].values)
    shp_reach_Q = np.zeros((rmax,))*np.nan
    lwidth      = np.zeros((rmax,))*np.nan
    Q_upper     = np.percentile(df_Q.values,95)
    Q_lower     = np.percentile(df_Q.values,5)
    Qmax        = np.max(df_Q,axis=0).reset_index()
    Qmax.columns = ['ID','Q']
    Q_max       = np.zeros((rmax,))*np.nan
    upper_limit = 5
    for i in range(0,len(shp_reach_Q)):
        ind = np.where((shp_reach['ID'].values[i]==Q['ID']))[0]
        if len(ind)>0: 
            ind            = int(ind)
            shp_reach_Q[i] = float(Q['Q'].values[ind])
            l_width        = (float(shp_reach_Q[i]) - Q_lower)/(Q_upper - Q_lower)
            lwidth[i]      = np.nanmin([np.nanmax([0, l_width  * upper_limit]),upper_limit]) + 0.5
            Q_max[i]       = float(Qmax['Q'].values[ind])
    shp_reach['Q'] = shp_reach_Q
    shp_reach['Qmax'] = Q_max
    shp_reach['lw']= lwidth
    
    return shp_reach, Q_upper, Q_lower


def Sobek_2_shp_WL(df_WL, shp_node):
    # Link Sobek data to shapefile: Waterlevel
    df_WL          = df_WL.drop(['Unnamed: 0','Location:'], axis=1)
    WL             = df_WL.iloc[0].reset_index() # timestep 0 -> change accordingly
    WL.columns     = ['ID','WL']
    
    nmax           = len(shp_node['ID'].values)
    shp_node_WL    = np.zeros((nmax, 1))*np.nan
    for i in range(0,len(shp_node_WL)):
        ind = np.where((shp_node['ID'].values[i]==WL['ID']))[0]
        if len(ind)>0: 
            ind               = int(ind)
            shp_node_WL[i]    = WL['WL'].values[ind]                
    shp_node['WL']     = shp_node_WL
       
    return shp_node


def point_to_poly(shp_node, shp_poly):
    CODE = 'CODE'
    
    # Add max(WLmax) to shp_poly
    shp_node_max = shp_node.groupby(CODE)['WL'].max().reset_index()
    WL_max_peil  = np.zeros(len(shp_poly['CODE']),)-99
    for i in range(0,len(shp_poly['CODE'])):
        idx            = np.where(shp_poly['CODE'][i]==shp_node_max[CODE])[0]
        if len(idx)>0:
            WL_max_peil[i] = shp_node_max.iloc[idx]['WL'].values
    shp_poly['WLmax']  = WL_max_peil
    
    return shp_poly

if __name__ == "__main__":
    
    # Input data
    shp_peil     = gpd.read_file(paths.shp_peil)
    shp_reach    = gpd.read_file(paths.shp_reach)
    shp_node     = gpd.read_file(paths.shp_node)
    df_Q         = pd.read_csv(paths.df_Q)
    df_WL        = pd.read_csv(paths.df_WL)
    
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
    shp_node.to_file(paths.Sobek_WL)          
    shp_reach.to_file(paths.Sobek_Q)          
    shp_peil.to_file(paths.shp_peil_WL)                
    
    