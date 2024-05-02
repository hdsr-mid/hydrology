# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:20:41 2023

@author: PetraH
"""


import os
import pandas as pd

    
if __name__ == "__main__":
    suffices = ['init','2u_GLG','24u_GLG','96u_GLG','2u_GHG','24u_GHG','96u_GHG']
    folder   = 'Sobek output' 
    for suffix in suffices:    
        print(suffix)
        df_Q      = pd.read_csv(folder + '/Q_'+suffix+'.csv',skiprows=4)
        df_WL     = pd.read_csv(folder + '/WL_'+suffix+'.csv',skiprows=4)
        
        if suffix == 'init':
            df_Q = df_Q[1:2]
            df_WL = df_WL[1:2]
        
        df_Q.to_csv(folder + '/Q_'+suffix+'.txt')
        df_WL.to_csv(folder + '/WL_'+suffix+'.txt')

    df_bui_GLG = pd.read_excel(folder + '/Bui.xlsx',sheet_name='GLG')
    df_bui_GHG = pd.read_excel(folder + '/Bui.xlsx',sheet_name='GHG')
    
    df_2u_GLG = df_bui_GLG[['P [mm/h]','2u']]
    df_24u_GLG = df_bui_GLG[['P [mm/h]','24u']]
    df_96u_GLG = df_bui_GLG[['P [mm/h]','96u']]
    df_2u_GLG.columns = ['Tijd','P [mm/h]']
    df_24u_GLG.columns = ['Tijd','P [mm/h]']
    df_96u_GLG.columns = ['Tijd','P [mm/h]']
    
    df_2u_GHG = df_bui_GHG[['P [mm/h]','2u']]
    df_24u_GHG = df_bui_GHG[['P [mm/h]','24u']]
    df_96u_GHG = df_bui_GHG[['P [mm/h]','96u']]
    df_2u_GHG.columns = ['Tijd','P [mm/h]']
    df_24u_GHG.columns = ['Tijd','P [mm/h]']
    df_96u_GHG.columns = ['Tijd','P [mm/h]']
    
    df_2u_GLG.to_csv(folder + '/Bui_init.txt')
    
    df_2u_GLG.to_csv(folder + '/Bui_2u_GLG.txt')
    df_24u_GLG.to_csv(folder + '/Bui_24u_GLG.txt')
    df_96u_GLG.to_csv('Bui_96u_GLG.txt')
    
    df_2u_GHG.to_csv(folder + '/Bui_2u_GHG.txt')
    df_24u_GHG.to_csv(folder + '/Bui_24u_GHG.txt')
    df_96u_GHG.to_csv(folder + '/Bui_96u_GHG.txt')
    