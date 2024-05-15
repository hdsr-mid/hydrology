# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:24:13 2023

@author: P. Hulsman
"""

# import sys
import os
import pandas as pd
import geopandas as gpd
# import rioxarray as rio
# import xarray as xr
import numpy as np
# import matplotlib.pyplot as plt
from pathlib import Path
import rasterio
from rasterstats import zonal_stats
from rasterio.fill import fillnodata
from rasterio.mask import mask
from scipy.ndimage import convolve
from shapely.geometry import mapping
import glob
from tqdm import tqdm

# Geef id van de zone op waarbinnen de percentages berekend moet worden
zonalid = 'OBJECTID'

# functie om percentielen van hoogtemodel te berekenen
def compute_zonal_stats(fp_vector, fp_raster, stats):
    """ Zonal stats voor AHN"""
    gdf = gpd.read_file(fp_vector)

    stats = zonal_stats(
        gdf,
        fp_raster,
        stats=stats
    )

    gdf_sel = gdf[[zonalid, 'geometry']]
    df_concat = pd.concat((gdf_sel, pd.DataFrame(stats)), axis=1)
    return df_concat

#    functie om uitsnede bodemhoogtemodel op te slaan
def get_updated_meta(f, out_img, out_transform):
    """ Raster meta-informatie"""
    out_meta = f.meta.copy()
    out_meta['transform'] = out_transform
    out_meta['compress'] = 'deflate'
    out_meta['height'] = out_img.shape[1]
    out_meta['width'] = out_img.shape[2]
    #out_meta['crs'] = {'init': 'epsg:28992'}
    out_meta['crs'] = None
    return out_meta

def fun_bodemberging(file_bodem_berging,GIS_path,GWST,w_init):
    '''bereken bodemberging
    classificatie:
            0: 0 - 5 cm     -> 0.05 m thick
            1: 5 - 15 cm    -> 0.10 m thick
            2: 15 - 30 cm   -> 0.15 m thick
            3: 30 - 60 cm   -> 0.30 m thick
            4: 60 - 100 cm  -> 0.40 m thick
            5: 100 - 200 cm -> 1.00 m thick
            aanname: > 200 cm heeft dezelfde porositeit als laag 5
            '''                
    
    # lees tif bestanden in
    SM_0 = rasterio.open(os.path.join(GIS_path,'HiHydro','HiHydro_wsat-'+w_init+'_0_clipped.tif'))
    SM_1 = rasterio.open(os.path.join(GIS_path,'HiHydro','HiHydro_wsat-'+w_init+'_1_clipped.tif'))
    SM_2 = rasterio.open(os.path.join(GIS_path,'HiHydro','HiHydro_wsat-'+w_init+'_2_clipped.tif'))
    SM_3 = rasterio.open(os.path.join(GIS_path,'HiHydro','HiHydro_wsat-'+w_init+'_3_clipped.tif'))
    SM_4 = rasterio.open(os.path.join(GIS_path,'HiHydro','HiHydro_wsat-'+w_init+'_4_clipped.tif'))
    SM_5 = rasterio.open(os.path.join(GIS_path,'HiHydro','HiHydro_wsat-'+w_init+'_5_clipped.tif'))  
    
    berging_l0_l0 = SM_0.read()*0.05
    berging_l0_l1 = SM_1.read()*0.10 + berging_l0_l0
    berging_l0_l2 = SM_2.read()*0.15 + berging_l0_l1
    berging_l0_l3 = SM_3.read()*0.30 + berging_l0_l2
    berging_l0_l4 = SM_4.read()*0.40 + berging_l0_l3
    berging_l0_l5 = SM_5.read()*1.00 + berging_l0_l4
    
    bodem_berging = np.zeros(GWST.read().shape) # output in [m]
    bodem_berging = np.where((GWST.read()>0)    & (GWST.read()<=0.05),                 SM_0.read()*GWST.read()         ,bodem_berging)
    bodem_berging = np.where((GWST.read()>0.05) & (GWST.read()<=0.15), berging_l0_l0 + SM_1.read()*(GWST.read() - 0.05),bodem_berging)
    bodem_berging = np.where((GWST.read()>0.15) & (GWST.read()<=0.3) , berging_l0_l1 + SM_2.read()*(GWST.read() - 0.15),bodem_berging)
    bodem_berging = np.where((GWST.read()>0.3)  & (GWST.read()<=0.6) , berging_l0_l2 + SM_3.read()*(GWST.read() - 0.30),bodem_berging)
    bodem_berging = np.where((GWST.read()>0.6)  & (GWST.read()<=1)   , berging_l0_l3 + SM_4.read()*(GWST.read() - 0.60),bodem_berging)        
    bodem_berging = np.where((GWST.read()>1)    & (GWST.read()<=2)   , berging_l0_l4 + SM_5.read()*(GWST.read() - 1.00),bodem_berging)
    bodem_berging = np.where((GWST.read()>2)                         , berging_l0_l5 + SM_5.read()*(GWST.read() - 2.00),bodem_berging)
    profile = GWST.profile
    with rasterio.open(file_bodem_berging, 'w', **profile) as dest:
            dest.write(bodem_berging)
    
def WaterBerging(GIS_path=None, temp_path=None, export_path=None,iteration=None, i_peilgebieden=None):
    
    # stats instellingen voor functie percentielen
    stats_water     = ['percentile_10','count']
    stats_verhard   = ['count']
    stats_onverhard = ['median','count','sum']

    # 3x3 window voor focal mean
    weights = np.ones((3,3))

    # Lees shp bestanden in geopandas
    bodem        = gpd.read_file(os.path.join(GIS_path,'Bodemkaart_HDSR.shp'))
    peilgebieden = gpd.read_file(os.path.join(GIS_path,'BR_Peilgebieden_20210819.shp'))
    insteek      = gpd.read_file(os.path.join(GIS_path,'BR_Insteekvlak.shp'))
    watervlak    = gpd.read_file(os.path.join(GIS_path,'BR_Watervlak.shp'))
            
    # lees tif bestanden in
    topo = rasterio.open(os.path.join(GIS_path,'AHN4 grondfilter.tif'))
    GWST = rasterio.open(os.path.join(GIS_path,'GxG','GWST_m.tif'))
    GHG  = rasterio.open(os.path.join(GIS_path,'GxG','GHG_m.tif'))
    GLG  = rasterio.open(os.path.join(GIS_path,'GxG','GLG_m.tif'))
    cellsize_topo = topo.transform[0] # m!
    cellsize_GW   = GWST.transform[0] # m!
    
    # bereken bodemberging
    file_bodem_berging_flc = os.path.join(GIS_path,'bodem_berging_wsat_wflc.tif')
    file_bodem_berging_wlp = os.path.join(GIS_path,'bodem_berging_wsat_wwlp.tif')
    file_bodem_berging_res = os.path.join(GIS_path,'bodem_berging_wsat_wres.tif')
    if not os.path.exists(file_bodem_berging_flc): fun_bodemberging(file_bodem_berging_flc,GIS_path,GHG,'wflc')
    if not os.path.exists(file_bodem_berging_wlp): fun_bodemberging(file_bodem_berging_wlp,GIS_path,GWST,'wwlp')
    if not os.path.exists(file_bodem_berging_res): fun_bodemberging(file_bodem_berging_res,GIS_path,GLG,'wres')
    
    # selecteer onverharde gebieden uit bodemeenheden
    # onverhard = bodem[(bodem['GRONDSOORT'] != 'Water') & (bodem['GRONDSOORT'] != 'Bebouwing')]
    onverhard = bodem[bodem['GRONDSOORT'] != 'Bebouwing']
    onverhard = onverhard.reset_index(drop=True)
    
    # selecteer verharde gebieden uit bodemeenheden
    verhard = bodem[bodem['GRONDSOORT'] == 'Bebouwing']
    verhard = verhard.reset_index(drop=True)

    # proces per peilgebied
    for i in tqdm(range(0,len(peilgebieden))):   
        # if i >= len(peilgebieden): continue
        
        idloc = int(peilgebieden.loc[i, zonalid])
        row = peilgebieden.iloc[[i]]
        
        # print('Iteratie: ', str(iteration), ', i: ', str(i),', Eenheid: ', str(idloc))
        # if idloc != 358: continue
        
        geometries = row.geometry.apply(mapping)

        # clip hoogtemodel op geselecteerde peilgebied
        cropped_ahn, out_transform = mask(topo, geometries, crop=True)
        out_meta = get_updated_meta(topo, cropped_ahn, out_transform)            
        file_cropped_ahn = os.path.join(temp_path, 'cropped_ahn_' + str(idloc) + '.tif')
        with rasterio.open(file_cropped_ahn, 'w', **out_meta) as dest:
            dest.write(cropped_ahn)

        # # clip GWST op geselecteerde peilgebied
        # cropped_GWST, out_transform_GWST = mask(GWST, geometries, crop=True)
        # out_meta_GWST = get_updated_meta(GWST, cropped_GWST, out_transform_GWST)
        # file_cropped_GWST = os.path.join(temp_path, 'cropped_GWST_' + str(idloc) + '.tif')
        # with rasterio.open(file_cropped_GWST, 'w', **out_meta_GWST) as dest:
        #     dest.write(cropped_GWST)            
        
        # Vull nodata op in geclipt hoogtemdel
        with rasterio.open(file_cropped_ahn) as src:
            profile = src.profile
            # arr = src.read(1)
            cropped_ahn = np.resize(cropped_ahn, (cropped_ahn.shape[1], cropped_ahn.shape[2]))
            arr_filled = fillnodata(cropped_ahn, mask=src.read_masks(1), max_search_distance=1000, smoothing_iterations=0)

        # In een window van 3x3 wordt een gemiddelde waarde uitgerekend op geclipte hoogtemodel
        ahnfilled = os.path.join(temp_path, 'ahnfilled_' + str(idloc) + '.tif')
        focal_mean = convolve(arr_filled, weights) / np.sum(weights)

        # resultaat focal mean opslaan
        with rasterio.open(ahnfilled, 'w', **profile) as dest:
            dest.write_band(1,focal_mean)
        
        # selecteer BGT-water in zone peilgebied
        out_differ_water = os.path.join(temp_path,'differ_water_' + str(idloc) + '.shp')
        gdf_differ_water = gpd.overlay(row, watervlak[['geometry']], how='intersection')
        gdf_differ_water.to_file(out_differ_water)
        
        # selecteer insteek in zone peilgebied
        out_differ_insteek = os.path.join(temp_path,'differ_insteek_' + str(idloc) + '.shp')
        gdf_differ_insteek = gpd.overlay(row, insteek[['geometry']], how='intersection')
        gdf_differ_insteek.to_file(out_differ_insteek)
        
        # selecteer BGT-onverhard in zone peilgebied
        out_differ_onverhard = os.path.join(temp_path,'differ_onverhard_' + str(idloc) + '.shp')
        gdf_differ_onverhard = gpd.overlay(row, onverhard[['geometry','BODEMTYPE','GRONDSOORT','SUBGROEP']], how='intersection')
        gdf_differ_onverhard.to_file(out_differ_onverhard)            
        
        # selecteer BGT-onverhard excl. insteek gebied in zone peilgebied
        out_differ_onverhard_excl_insteek = os.path.join(temp_path,'differ_onverhard_excl_insteek_' + str(idloc) + '.shp')                
        if len(gdf_differ_onverhard)>0:
            gdf_differ_onverhard_excl_insteek = gpd.overlay(gdf_differ_onverhard, insteek[['geometry']], how='difference')
            gdf_differ_onverhard_excl_insteek.to_file(out_differ_onverhard_excl_insteek)            
        
        # selecteer BGT-onverhard in zone peilgebied
        out_differ_verhard = os.path.join(temp_path,'differ_verhard_' + str(idloc) + '.shp')
        gdf_differ_verhard = gpd.overlay(row, verhard[['geometry']], how='intersection')
        gdf_differ_verhard.to_file(out_differ_verhard)

        # referentie peil
        PEIL              = np.where(row["ZOMERPEIL"]==0,row["VASTPEIL"],row["ZOMERPEIL"])
        PEIL              = np.where(PEIL==-999,np.nan,PEIL)
        
        '''Berging oppervlakte water'''
        if (len(gdf_differ_water)>0) & (len(gdf_differ_onverhard)>0): 
            # Verkrijg statistiek binnen de gedefinieerde zones
            gdf_zonal_stats = compute_zonal_stats(out_differ_onverhard_excl_insteek, ahnfilled, stats_water)

            # Verwijder, voeg kolommen toe en bereken per zone een maaiveldcurve: hoogte(m), oppervlak(m2) en volume(m3)
            df1 = pd.DataFrame(gdf_zonal_stats.drop(columns='geometry'))
            df1["hoogte_10p"] = df1["percentile_10"]
            df1["A_insteek"]  = gdf_differ_insteek['geometry'].area.astype(float)
            df1["A_water"]    = gdf_differ_water['geometry'].area.astype(float)
            df1["PEIL"]       = float(PEIL)
            df1["ruimte_ver"] = df1["hoogte_10p"] - df1["PEIL"]
            df1["ruimte_ver"] = np.where(df1["ruimte_ver"]<0,0,df1["ruimte_ver"])
            # df1["ruimte_ver"] = np.where(df1["ruimte_ver"]<0.1,0.1,df1["ruimte_ver"])
            df1["volume"]     = df1["ruimte_ver"]*(0.5 * (df1["A_insteek"] + df1["A_water"]) )
            df_water          = df1.copy()                     
            del df1,gdf_zonal_stats                       
        else:
            df_water           = row[[row.columns[0]]]
            df_water["volume"] = np.nan    
            
        '''Berging verhard gebied (bebouwd)'''
        if len(gdf_differ_verhard)>0: 
            # Verkrijg statistiek binnen de gedefinieerde zones
            gdf_zonal_stats = compute_zonal_stats(out_differ_verhard, ahnfilled, stats_verhard)

            # Verwijder, voeg kolommen toe en bereken per zone een maaiveldcurve: hoogte(m), oppervlak(m2) en volume(m3)
            df1 = pd.DataFrame(gdf_zonal_stats.drop(columns='geometry'))
            df1["A_verhard"]  = gdf_differ_verhard['geometry'].area.astype(float)                
            df1["ruimte_ver"] = 5 / 1000 # aanname gemiddelde berging riool [m]
            df1["volume"]     = df1["ruimte_ver"]*df1["A_verhard"]                
            df_verhard        = df1.copy()                
            del df1,gdf_zonal_stats
        else:
            df_verhard           = row[[row.columns[0]]]
            df_verhard["volume"] = np.nan
            
        '''Berging onverhard gebied (bodem)'''            
        if len(gdf_differ_onverhard)>0:                
            # Verkrijg statistiek binnen de gedefinieerde zones
            gdf_zonal_stats_flc = compute_zonal_stats(out_differ_onverhard, file_bodem_berging_flc, stats_onverhard)
            gdf_zonal_stats_wlp = compute_zonal_stats(out_differ_onverhard, file_bodem_berging_wlp, stats_onverhard)
            gdf_zonal_stats_res = compute_zonal_stats(out_differ_onverhard, file_bodem_berging_res, stats_onverhard)
            
            # Verwijder, voeg kolommen toe en bereken per zone een maaiveldcurve: hoogte(m), oppervlak(m2) en volume(m3)
            df1 = pd.DataFrame(gdf_zonal_stats_flc.drop(columns='geometry'))  
            df1["A_onverhard"] = gdf_differ_onverhard['geometry'].area.astype(float)
            df1['volume']      = df1["A_onverhard"]*df1['median'] 
            df_onverhard_flc   = df1.copy() 
            del df1,gdf_zonal_stats_flc                     
            
            df1 = pd.DataFrame(gdf_zonal_stats_wlp.drop(columns='geometry'))  
            df1["A_onverhard"] = gdf_differ_onverhard['geometry'].area.astype(float)
            df1['volume']      = df1["A_onverhard"]*df1['median'] 
            df_onverhard_wlp   = df1.copy() 
            del df1,gdf_zonal_stats_wlp                     
            
            df1 = pd.DataFrame(gdf_zonal_stats_res.drop(columns='geometry'))  
            df1["A_onverhard"] = gdf_differ_onverhard['geometry'].area.astype(float)
            df1['volume']      = df1["A_onverhard"]*df1['median'] 
            df_onverhard_res   = df1.copy() 
            del df1,gdf_zonal_stats_res                     
        else:
            df_onverhard_flc           = row[[row.columns[0]]]
            df_onverhard_wlp           = row[[row.columns[0]]]
            df_onverhard_res           = row[[row.columns[0]]]
            df_onverhard_flc["volume"] = np.nan    
            df_onverhard_wlp["volume"] = np.nan    
            df_onverhard_res["volume"] = np.nan    
                    
        '''Gebieden samenvoegen'''
        df1               = row[[row.columns[0]]]
        df1["area"]       = float(row['geometry'].area)
        df1["BWm3"]     = np.nansum(df_water['volume'])     # som van alle water oppervlaktes binnen 1 peilgebied
        df1["BVm3"]     = np.nansum(df_verhard['volume'])
        df1["BOm3_min"]= np.nansum(df_onverhard_flc['volume'])
        df1["BOm3_gem"]= np.nansum(df_onverhard_wlp['volume'])
        df1["BOm3_max"]= np.nansum(df_onverhard_res['volume'])
        df1["Btotm3_min"]   = np.nansum(np.array([df1["BWm3"] + df1["BVm3"] + df1["BOm3_min"]]))
        df1["Btotm3_gem"]   = np.nansum(np.array([df1["BWm3"] + df1["BVm3"] + df1["BOm3_gem"]]))
        df1["Btotm3_max"]   = np.nansum(np.array([df1["BWm3"] + df1["BVm3"] + df1["BOm3_max"]]))
        
        df1["BWmm"]     = np.nansum(df_water['volume'])  / df1["area"] * 1000 # m -> mm
        df1["BVmm"]     = np.nansum(df_verhard['volume'])  / df1["area"] * 1000 # m -> mm
        df1["BOmm_min"]= np.nansum(df_onverhard_flc['volume'])  / df1["area"] * 1000 # m -> mm
        df1["BOmm_gem"]= np.nansum(df_onverhard_wlp['volume'])  / df1["area"] * 1000 # m -> mm
        df1["BOmm_max"]= np.nansum(df_onverhard_res['volume'])  / df1["area"] * 1000 # m -> mm
        df1["Btotmm_min"]   = df1["Btotm3_min"] / df1["area"] * 1000 # m -> mm
        df1["Btotmm_gem"]   = df1["Btotm3_gem"] / df1["area"] * 1000 # m -> mm
        df1["Btotmm_max"]   = df1["Btotm3_max"] / df1["area"] * 1000 # m -> mm
        
        # Resultaat samenvoegen met oorsponkelijke shape en opslaan als shapefile
        resultdef            = pd.merge(peilgebieden, df1, on=[zonalid, zonalid])
        output_file = os.path.join(temp_path,"maaiveldcurve_" + str(idloc) + ".shp")
        resultdef.to_file(output_file)

        # opschonen data
        os.remove(file_cropped_ahn)
        os.remove(ahnfilled)
        # os.remove(file_cropped_GWST)
        [os.remove(os.path.join(temp_path,f)) for f in [file for file in os.listdir(temp_path) if file.startswith('differ') and str(idloc) in file]]            

def fun_combine_shp(GIS_path=None, temp_path=None, export_path=None):
    # haal maaiveldcurve shapes op
    files = glob.glob(temp_path+'/maaiveldcurve_*.shp')

    # maak een lege lijst om dataframes in op te slaan
    li = []

    # loop through list of files and read each one into a dataframe and append to list
    for f in files:
        temp_df = gpd.read_file(f)
        li.append(temp_df)

    # concate alle dataframes
    df = pd.concat(li, axis=0)
    columns_final = ['OBJECTID','area', 'geometry',
                      'BWm3','BVm3','BOm3_min','BOm3_gem','BOm3_max','Btotm3_min','Btotm3_gem','Btotm3_max',
                      'BWmm','BVmm','BOmm_min','BOmm_gem','BOmm_max','Btotmm_min','Btotmm_gem','Btotmm_max']
    df.drop([col for col in df.columns if col not in columns_final], axis=1, inplace=True)        
    [os.remove(os.path.join(temp_path,f)) for f in [file for file in os.listdir(temp_path) if file.startswith('maaiveldcurve')]]
    df = df.fillna(-999)
    
    # schrijf resultaat als shapefile
    df.to_file(os.path.join(export_path, 'hdsr_maaiveldcurves4.shp'))
        
def fun_afronden(GIS_path=None, temp_path=None, export_path=None):
    file = os.path.join(export_path, 'hdsr_maaiveldcurves4.shp')
    df   = gpd.read_file(file)
    
    columns_final = ['area', 'BWm3','BVm3','BOm3_min','BOm3_gem','BOm3_max','Btotm3_min','Btotm3_gem','Btotm3_max',
                      'BWmm','BVmm','BOmm_min','BOmm_gem','BOmm_max','Btotmm_min','Btotmm_gem','Btotmm_max']
    
    for col in columns_final:
        df[col]   = np.round(df[col])
    
    df.to_file(os.path.join(export_path, 'bergingscapaciteit.shp'))
    
