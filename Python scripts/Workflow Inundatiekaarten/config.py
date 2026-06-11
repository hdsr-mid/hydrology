# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 14:22:39 2023

@author: PetraH
"""

import os
import sys

# define variable names
CODE = 'CODE' # column name for CODE in peilgebieden shapefile (it seems to change...)

# define all paths here
root   = os.getcwd()

# Paths: Input files
raster_file = os.path.join(root,'GIS_input', "AHN3_ruw_clip.tif")
shp_node    = os.path.join(root,"GIS_input",'3b_nod_edited_clip.shp')
shp_peil    = os.path.join(root,'GIS_input','BR-VS_Peilgebieden_clip.shp')
shp_extent  = os.path.join(root,'GIS_input','EvS_gebied.shp')
shp_insteek = os.path.join(root,'GIS_input','BR_VS_Insteekvlak.shp')
shp_water   = os.path.join(root,'GIS_input','BR_VS_Watervlak.shp')
shp_afvoer  = os.path.join(root,'GIS_input','BR_VS_Afvoergebieden_clip.shp')
shp_afw     = os.path.join(root,'GIS_input','Afwateringseenheden_2020_clip.shp')
df_WL       = os.path.join(root,'Sobek output','WL.txt')

# Paths: output files
output      = os.path.join(root,'output')
GIS_output  = os.path.join(root,'GIS_output')
shp_Sobek_WL= os.path.join(GIS_output,'Sobek_WLmax.shp')
shp_peil_WLmax = os.path.join(GIS_output,'BR-VS_Peilgebieden_WLmax.shp')
raster_inun    = os.path.join(GIS_output,"inundatiekaart_cm.tif")
gpk_out        = os.path.join(root, 'GIS_output','inundatie.gpkg') 

# Paths: temp files
GIS_temp    = os.path.join(root,'GIS_temp')
shp_poly_10p= os.path.join(GIS_temp,'BR-VS_Peilgebieden_clip_edited_10p.shp')
raster_peil = os.path.join(GIS_temp,'Peil.tif')
raster_clip = os.path.join(GIS_temp,"AHN3_ruw_clip.tif")
raster_prep = os.path.join(GIS_temp,"AHN3_ruw_clip_prep.tif")
shp_droognat  = os.path.join(GIS_temp,'shp_droognat.shp')    
tif_clip      = os.path.join(GIS_temp,'raster_clip.tif')
tif_clip_binar= os.path.join(GIS_temp,'raster_clip_binary.tif')
shp_out_inun  = os.path.join(GIS_temp,'inundatie.shp')    
shp_out_nat   = os.path.join(GIS_temp,'nat.shp')    
    
    