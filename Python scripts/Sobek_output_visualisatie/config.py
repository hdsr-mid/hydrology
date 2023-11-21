# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 14:22:39 2023

@author: PetraH
"""

import os
import sys

# define all paths here

root   = os.getcwd()

# Paths: Input files
raster_file = os.path.join(root,'GIS_input', "AHN3_ruw.tif")
shp_peil    = os.path.join(root,'GIS_input','BR-VS_Peilgebieden_clip_edited.shp')
# BR_peil     = os.path.join(root,'GIS_input','BR_VS_Peilgebieden.shp')
shp_reach   = os.path.join(root,"GIS_input","RchSegments.shp")
shp_node    = os.path.join(root,"GIS_input",'3b_nod_edited_joined_drooglegging.shp')
shp_extent  = os.path.join(root,'GIS_input','EvS_gebied_dissolved.shp')
shp_water   = os.path.join(root,'GIS_input','BR_VS_Watervlak.shp')
shp_insteek = os.path.join(root,'GIS_input','BR_VS_Insteekvlak.shp')
df_Q        = os.path.join(root,'Sobek output','Q.txt')
df_WL       = os.path.join(root,'Sobek output','WL.txt')
df_P        = os.path.join(root,'Sobek output','Bui.txt')

# Paths: output files
output      = os.path.join(root,'output')
GIS_output  = os.path.join(root,'GIS_output')
GIS_temp    = os.path.join(root,'GIS_temp')
Sobek_dWL   =os.path.join(GIS_output,'Sobek_dWL.shp')
peil_clip   = os.path.join(GIS_temp,'BR_VS_Peilgebieden_clip.shp')
raster_peil = os.path.join(GIS_temp,'Peil.tif')
raster_clip = os.path.join(GIS_temp,"AHN3_ruw_clip.tif")
raster_prep = os.path.join(GIS_temp,"AHN3_ruw_clip_prep.tif")
raster_inun = os.path.join(GIS_output,"inundatiekaart.tif")

