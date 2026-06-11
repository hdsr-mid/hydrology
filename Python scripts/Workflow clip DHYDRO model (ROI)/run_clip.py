# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21  2025

@author: MaaS-user

Stappen:
    1. clip netcdf
    2. clip crsloc.ini o.b.v. netcdf
    3. clip crsdef.ini o.b.v. crsloc.ini
    4. clip overige ini bestanden o.b.v. netcdf
    
"""

import os
from fun_clip_init import clip_ini
from fun_clip_netcdf import clip_netcdf
import shutil
from fun_round import round_netcdf, round_pliz

if __name__ == "__main__":
    root  = r'Z:\Petra\DHydro'
    model = 'Combined_V3.16_100m_'
    
    mtype     = '1D2D'
    scenarios = ['_1D2Ddetail','_1D2Driv','_1D2Drivplus']
    # scenarios = ['_1D2D_nolinks_detail','_1D2D_nolinks_riv']
    # mtype     = '1D'
    # scenarios = ['_1Ddetail','_1Driv']
    
    # afronden, nieuw basismodel
    s = os.path.join(root, model + mtype + '_newsettings')  # source directory
    d = os.path.join(root, model + mtype + '_rounded')  # destination directory
    if not os.path.exists(d):
        shutil.copytree(s, d)
        round_netcdf(os.path.join(d,'dflowfm'))
        # round_pliz(os.path.join(d,'dflowfm'))
    else:
        print("Already exists")
    
    
    for scenario in scenarios:
        print(scenario)
        
        # create new folder
        s = os.path.join(root, model + mtype + '_rounded')  # source directory
        d = os.path.join(root, model + mtype + '_rounded' + scenario)  # destination directory
        if not os.path.exists(d):
            shutil.copytree(s, d)
        else:
            print("Already exists")
        
        path  = os.path.join(root,d,'dflowfm')
        clip_netcdf(path)
        
        clip_ini(path,'structures.ini','[Structure]\n')
        clip_ini(path,'crsloc.ini','[CrossSection]\n')
        clip_ini(path,'crsdef.ini', '[Definition]\n')
        clip_ini(path,'InitialWaterLevel.ini','[Branch]\n')    
        clip_ini(path,'roughness-FloodPlain1.ini','[Branch]\n')
        clip_ini(path,'roughness-Main.ini','[Branch]\n')
    
    
    
