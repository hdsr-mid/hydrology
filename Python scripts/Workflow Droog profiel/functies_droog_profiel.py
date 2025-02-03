# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:48:37 2024

@author: PetraH
"""

import numpy as np
import geopandas as gpd
import pandas as pd
from osgeo import gdal 
from shapely.geometry import Point
from shapely.geometry import LineString 
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")

AHN_noData = -1E+22
d_extrapolate = 15 # extrapoleer afstand voorbij uiterst profielpunt; deze moet minimaal voorbij de insteeklijn gaan

def format_WIT(paths):
    # Input data
    gdf = gpd.read_file(paths.shp_WIT)    
    
    # format shapefile
    gdf['punttype']   = [p.replace('_','') for p in gdf['punt_typ_i'].values]
    gdf['profielmet'] = gdf['pm_id'].astype(str)    
    gdf['datum']      = [str(s).split('_')[-1] for s in gdf['pl_ids'].values]
    gdf               = gdf[['profielmet','datum','punttype','afstand','slibhoogte','X','Y','geometry']]
    
    # add Hydro Object code
    gdf = add_HO_ID(paths, gdf)
    
    # keep only profiles in primary waterways
    names_all= np.unique(gdf['profielmet'].values)
    gdf      = gdf[gdf['HO_cat']==1]    
    names_sel= np.unique(gdf['profielmet'].values)
    names   = [n for n in names_all if n not in names_sel]
    with open(paths.txt_dropped, 'w') as outfile:
        outfile.write('Profiles dropped because they are not located in primary waterways:\n')
        outfile.write(str(names))
        outfile.write('\n')    
    print(len(names), ' of ',len(names_all), ' profiles are dropped because they are not located in primary waterways. Check ', paths.txt_dropped, ' for the ID codes.')
    
    # keep only profiles close to a waterway
    names_all= np.unique(gdf['profielmet'].values)
    gdf      = gdf[gdf['HO_afstd']<1]    
    names_sel= np.unique(gdf['profielmet'].values)
    names   = [n for n in names_all if n not in names_sel]
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because their min distance from the waterway is more than 1 m:\n')
        outfile.write(str(names))
        outfile.write('\n')    
    print(len(names), ' of ',len(names_all), ' profiles are dropped because their min distance from the waterway is more than 1 m. Check ', paths.txt_dropped, ' for the ID codes.')
    
    # keep only profiles with identical codes for the hydroobject and insteekvlak
    names_all= np.unique(gdf['profielmet'].values)
    gdf      = gdf[gdf['HO_ID']==gdf['INSTK_ID']]    
    names_sel= np.unique(gdf['profielmet'].values)
    names   = [n for n in names_all if n not in names_sel]
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because their codes do not match (hydroobject vs. insteekvlak):\n')
        outfile.write(str(names))
        outfile.write('\n')    
    print(len(names), ' of ',len(names_all), ' profiles are dropped because their codes do not match (hydroobject vs. insteekvlak). Check ', paths.txt_dropped, ' for the ID codes.')
    
    # Save
    gdf.to_file(paths.shp_WIT_edited)


def add_HO_ID(paths, gdf):
    ''' Add the code and smallest distance to a Hydro object to the profiles'''
    
    # Input data
    shp_HO    = gpd.read_file(paths.shp_HO)
    shp_HO    = shp_HO[['CODE','CATEGORIEO','geometry']]
    shp_INSTK = gpd.read_file(paths.shp_insteek)
    shp_INSTK = shp_INSTK[['CODE','CODE_HO','geometry']]
    
    # Group points based on the column 'profielmet'
    grouped = gdf.groupby('profielmet')
    
    # Loop through each profielmet-code
    for name, group in grouped:    
        ### Hydro object ###
        # Spatial join
        gdf_joined = gpd.sjoin_nearest(group, shp_HO, distance_col="distances")
        
        # get smallest distance
        distance = np.min(gdf_joined.distances.values)
        
        # add info to gdf        
        gdf.loc[(gdf['profielmet'] == name), 'HO_ID']      = gdf_joined.CODE.values[0]
        gdf.loc[(gdf['profielmet'] == name), 'HO_cat']     = gdf_joined.CATEGORIEO.values[0]
        gdf.loc[(gdf['profielmet'] == name), 'HO_afstd']   = distance
        
        ### Insteekvlak ###
        # Spatial join
        gdf_joined = gpd.sjoin_nearest(group, shp_INSTK, distance_col="distances")
        
        # get smallest distance
        distance = np.min(gdf_joined.distances.values)
        
        # add info to gdf        
        gdf.loc[(gdf['profielmet'] == name), 'INSTK_ID']      = gdf_joined.CODE_HO.values[0]
        gdf.loc[(gdf['profielmet'] == name), 'INSTK_afstd']   = distance
        
    return gdf  
    
def fun_preprocessing(paths):
    '''Remove profiles based on the coordinates for A22L and A22R'''
    
    # Input data
    gdf = gpd.read_file(paths.shp_WIT_edited)
    
    # Drop profiles if they don't have A22L AND A22R
    names    = np.unique(gdf['profielmet'].values)
    gdf_A22L = gdf[gdf['punttype']=='A22L']    
    gdf_A22R = gdf[gdf['punttype']=='A22R']    
    names_L  = gdf_A22L['profielmet'].values
    names_R  = gdf_A22R['profielmet'].values
    
    names   = [n for n in names if (n not in names_L) or (n not in names_R)]
    for name in names:
        gdf     = gdf.drop(gdf[gdf['profielmet']==name].index)
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because the code A22L or A22R is missing:\n')
        outfile.write(str(names))
        outfile.write('\n')
    print(len(names), ' of ',len(gdf_A22L), ' profiles are dropped because the code A22L or A22R is missing. Check ', paths.txt_dropped, ' for the ID codes.')
    
    # Drop profiles if the coordinates for A22L and A22R are the same
    gdf_A22L = gdf[gdf['punttype']=='A22L']    
    gdf_A22R = gdf[gdf['punttype']=='A22R']    
    diff_X = abs(np.round(gdf_A22L['X'].values-gdf_A22R['X'].values,2))
    diff_Y = abs(np.round(gdf_A22L['Y'].values-gdf_A22R['Y'].values,2))
    
    ind = np.where((diff_X + diff_Y == 0))[0]
    names = gdf_A22L.iloc[ind]['profielmet']
    for name in names:
        gdf     = gdf.drop(gdf[gdf['profielmet']==name].index)
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because the codes A22L and A22R have the same coordinates:\n')
        outfile.write(str(names))
        outfile.write('\n')
    print(len(names), ' of ',len(gdf_A22L), ' profiles are dropped because the codes A22L and A22R have the same coordinates. Check ', paths.txt_dropped, ' for the ID codes.')
    
    # save file
    gdf.to_file(paths.shp_profielen_nat)
    
def fun_extrapolate(paths):
    ''' Add points at the edges of each profile covering the dry part of the profile
    These points are temporary and slightly beyond the profile edge
    These points are called A0 en Z0 to ensure they end up on the most left (A0) or right (Z0) of the profile considering other profile codes'''
    
    # Input data
    gdf = gpd.read_file(paths.shp_profielen_nat)
    df  = gpd.GeoDataFrame(columns = ['profielmet','punttype','HO_ID', 'HO_afstd','HO_cat','X', 'Y'])
    
    # Iterate through each profile
    for name, group in tqdm(gdf.groupby('profielmet')): 
        # select profile edge & interpolate
        group_edge = group[(group['punttype']=='A22L') | (group['punttype']=='A22R')]    
        x = group_edge['X']
        y = group_edge['Y']
        
        # regression line
        fy = np.poly1d(np.polyfit(x, y, 1)) # y = ax + b
        fx = np.poly1d(np.polyfit(y, x, 1)) # x = ay + b
        
        # extrapolation coordinates (with a maximum distance of d_extrapolate)
        dx = np.diff(group_edge['X'].values)
        if dx > 0: dx = d_extrapolate 
        if dx < 0: dx = -d_extrapolate 
        dy = np.diff(group_edge['Y'].values)
        if dy > 0: dy = d_extrapolate
        if dy < 0: dy = -d_extrapolate
        
        # Left side
        punttype = 'A0'
        if abs(fy[1]) < 1.5: # not a vertical line (based on the slope of the regression line)
            # estimate y-coord from regression (assuming the line is not vertical)
            Xcoord   = group_edge[group_edge['punttype']=='A22L']['X'].values - dx
            Ycoord   = fy(Xcoord)
        else:
            # estimate x-coord from regression (assuming line is vertical)
            Ycoord   = group_edge[group_edge['punttype']=='A22L']['Y'].values - dy
            Xcoord   = fx(Ycoord)
        # add point to dataframe
        df_new = pd.DataFrame({'profielmet': name,'punttype': punttype, 
                               'HO_ID': group['HO_ID'].values[0], 'HO_afstd': group['HO_afstd'].values[0],'HO_cat':group['HO_cat'].values[0], 
                               'X': Xcoord, 'Y': Ycoord})
        df = pd.concat([df, df_new])        
        
        # Right side
        punttype = 'Z0'        
        if abs(fy[1]) < 1.5: # not a vertical line (based on the slope of the regression line)
            # estimate y-coord from regression (assuming the line is not vertical)
            Xcoord   = group_edge[group_edge['punttype']=='A22R']['X'].values + dx
            Ycoord   = fy(Xcoord)            
        else:
            # estimate x-coord from regression (assuming line is vertical)
            Ycoord   = group_edge[group_edge['punttype']=='A22R']['Y'].values + dy
            Xcoord   = fx(Ycoord)
        df_new = pd.DataFrame({'profielmet': name,'punttype': punttype, 
                               'HO_ID': group['HO_ID'].values[0], 'HO_afstd': group['HO_afstd'].values[0],'HO_cat':group['HO_cat'].values[0],
                               'X': Xcoord, 'Y': Ycoord})
        df = pd.concat([df, df_new])   
        
    # save dataframe as shapefile
    gdf  = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.X,df.Y)])
    gdf.crs='EPSG:28992'
    gdf.to_file(paths.shp_points)
    
    
def points_2_line(paths):
    ''' Draw a line through the 2 extreme edges (A0 and Z0) of the profile''' 
    
    # Input data
    gdf = gpd.read_file(paths.shp_points)
    
    # Filter the points with code 'A0' and 'Z0' 
    selected_points = gdf[gdf['punttype'].isin(['A0', 'Z0'])] 
    
    # Iterate through profiles
    lines = [] 
    for profielmet, group in tqdm(selected_points.groupby('profielmet')): 
        if len(group) == 2: 
            point_A0 = group[group['punttype'] == 'A0'].geometry.iloc[0] 
            point_Z0 = group[group['punttype'] == 'Z0'].geometry.iloc[0] 
           
            # Create a line between these points 
            line = LineString([point_A0, point_Z0]) 
            
            # Append the line and its 'profielmet' code to the list 
            lines.append({'geometry': line, 'profielmet': profielmet,
                          'HO_ID': group['HO_ID'].values[0], 'HO_afstd': group['HO_afstd'].values[0],'HO_cat':group['HO_cat'].values[0]}) 
            
        else: 
            # Handle cases where a pair is not found 
            print(f"Pair not found for 'profielmet' code: {profielmet}") 
    
    # Create a GeoDataFrame with the lines
    lines_gdf = gpd.GeoDataFrame(lines) 
    # Set the CRS to match the input points 
    lines_gdf.crs = selected_points.crs 
    # Save to file 
    lines_gdf.to_file(paths.shp_lines) 
       
def clip_line_by_polygon(paths):
    ''' Clip line by polygon'''
    
    # Input data
    gdf_line    = gpd.read_file(paths.shp_lines)
    gdf_polygon = gpd.read_file(paths.shp_insteek)
    shp_HO      = gpd.read_file(paths.shp_HO)[['OBJECTID','CODE', 'geometry']]
    
    # Clip
    intersection = gpd.clip(gdf_line,gdf_polygon)
    
    # Remove multilines
    gdf_single        = intersection[intersection.geometry.geom_type == 'LineString']
    gdf_multi         = intersection[intersection.geometry.geom_type == 'MultiLineString'].explode()
    gdf_multi_joined  = gpd.sjoin(gdf_multi,shp_HO, how = 'left')
    gdf_multi_joined  = gdf_multi_joined[~gdf_multi_joined['index_right'].isna()]
    gdf_multi_joined  = gdf_multi_joined[gdf_multi_joined['HO_ID']==gdf_multi_joined['CODE']]
    
    # Remove profiels still consisting of multiple lines
    names             = gdf_multi_joined['profielmet'].drop_duplicates()
    gdf_multi_joined  = gdf_multi_joined.drop_duplicates(subset=['profielmet'],keep=False)
    names             = [n for n in names if n not in gdf_multi_joined['profielmet'].values]
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because they consist of multiple lines (e.g. at a bend):\n')
        outfile.write(str(names))
        outfile.write('\n')
    print(len(names), ' of ',len(names), ' profiles are dropped because they consist of multiple lines (e.g. at a bend). Check ', paths.txt_dropped, ' for the ID codes.')
    
    # merge remaining profiles
    gdf               = gpd.GeoDataFrame(pd.concat([gdf_single, gdf_multi_joined], ignore_index=True))
    
    # Save to file
    gdf.to_file(paths.shp_lines_insteek) 
    
def lines_2_points(paths):
    ''' Extract points at the edges of the line'''
    
    # Input data
    gdf_line    = gpd.read_file(paths.shp_lines_insteek)
      
    # Create new gdf
    points = []
    for i in tqdm(range(0,len(gdf_line))):
        line    = gdf_line.iloc[i]
        X_coord = line['geometry'].xy[0]
        Y_coord = line['geometry'].xy[1]
        for x in range(0,len(X_coord)):
            geom    = Point(X_coord[x],Y_coord[x])
            points.append({'profielmet': gdf_line['profielmet'].values[i], 'geometry': geom,
                           'HO_ID': gdf_line['HO_ID'].values[i], 'HO_afstd': gdf_line['HO_afstd'].values[i],'HO_cat':gdf_line['HO_cat'].values[i]}) 
    gdf_points = gpd.GeoDataFrame(points) 
    
    # Save shp file
    gdf_points.to_file(paths.shp_points_insteek) 
    
def get_elevation(x, y, gt, raster_band): 
    '''  Function to get the value from the raster ''' 
    px = int((x - gt[0]) / gt[1])  # x pixel 
    py = int((y - gt[3]) / gt[5])  # y pixel 
    try: 
        topo     = raster_band.ReadAsArray(px, py, 1, 1)[0][0] 
        if topo < AHN_noData:
            topo = [raster_band.ReadAsArray(px + dx, py + dy, 1, 1)[0][0] for dx in np.arange(-0.5,1,0.5) for dy in np.arange(-0.5,1,0.5)]
            topo = [t for t in topo if t > -1e22]
            if len(topo)==0: # still no data
                topo = AHN_noData
            elif np.array(topo).max() - np.array(topo).min() > 0.05: # too high height-differences to take the mean
                topo = AHN_noData
            else:
                topo = np.mean(np.array(topo))
        return topo
    except Exception as e: 
        print(f"Error reading elevation at point ({x}, {y}): {e}") 
      
   
def extract_data_along_points(paths):
    ''' Add raster data to polypoints''' 
    
    # Input data
    points = gpd.read_file(paths.shp_points_insteek) 
    
    # get AHN
    raster_path    = paths.tif_AHN
    raster         = gdal.Open(raster_path) 
    gt             = raster.GetGeoTransform() 
    raster_band    = raster.GetRasterBand(1) 
    points['AHN']  = points.apply(lambda row: get_elevation(row.geometry.x, row.geometry.y, gt, raster_band), axis=1) 
    
    # Save the GeoDataFrame with the elevation data to a new shapefile (optional)     
    points.to_file(paths.shp_points_insteek_AHN) 

    # Close the raster file 
    raster = None 
  
def merge_shapefiles(paths):
    ''' Merge shapefiles: add points of the dry part of the profile to the shapefile with the points in the wet part
    Update the distance column ('afstand') and point type ('punttype') column for the new points'''
    
    # Input data
    gdf1 = gpd.read_file(paths.shp_points_insteek_AHN)
    gdf2 = gpd.read_file(paths.shp_profielen_nat)
    gdf1['profielmet'] = gdf1['profielmet'].values.astype(str) # make sure it is a string
    gdf2['profielmet'] = gdf2['profielmet'].values.astype(str) # make sure it is a string
    code_links = 'A24'
    code_rechts= 'A28' 
    
    # Update de X and Y columns based on the existing geometry
    gdf1['X'] = gdf1.geometry.x
    gdf1['Y'] = gdf1.geometry.y
    gdf2['X'] = gdf2.geometry.x
    gdf2['Y'] = gdf2.geometry.y
    
    # Sometimes profiles are not located in the "insteekvlak" and are therefore not included in shp_points_insteek_AHN anymore
    # These profiles are now removed from shp_profielen_nat
    names_dry = np.unique(gdf1['profielmet'].values)
    names_wet = np.unique(gdf2['profielmet'].values)
    names_drop= [n for n in names_wet if n not in names_dry]
    for name in names_drop:
        gdf2     = gdf2.drop(gdf2[gdf2['profielmet']==name].index)
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because it is not located in the insteekvlak:\n')
        outfile.write(str(names_drop))
        outfile.write('\n')
    print(len(names_drop), ' of ',len(names_wet), ' profiles are dropped because it is not located in the insteekvlak. Check ', paths.txt_dropped, ' for the ID codes.')
    
    # Sometimes the WIT profiles are shifted compared to satellite imagery (and BR)
    # These profiles are removed
    X_insteek = gdf1.groupby('profielmet').agg({'X': ['min', 'max']}).reset_index()
    Y_insteek = gdf1.groupby('profielmet').agg({'Y': ['min', 'max']}).reset_index()
    X_nat     = gdf2.groupby('profielmet').agg({'X': ['min', 'max']}).reset_index()
    Y_nat     = gdf2.groupby('profielmet').agg({'Y': ['min', 'max']}).reset_index()
    ind_drop = np.where((np.round(X_insteek['X']['max'].values,2) < np.round(X_nat['X']['max'].values,2)) | 
             (np.round(X_insteek['X']['min'].values,2) > np.round(X_nat['X']['min'].values,2)) |
             (np.round(Y_insteek['Y']['max'].values,2) < np.round(Y_nat['Y']['max'].values,2)) |
             (np.round(Y_insteek['Y']['min'].values,2) > np.round(Y_nat['Y']['min'].values,2))
             )
    names = list(X_insteek.iloc[ind_drop]['profielmet'].values)
    for name in names:
        gdf1     = gdf1.drop(gdf1[gdf1['profielmet']==name].index)
        gdf2     = gdf2.drop(gdf2[gdf2['profielmet']==name].index)
    with open(paths.txt_dropped, 'a') as outfile:
        outfile.write('Profiles dropped because the wet profile is shifted such that A24/A28 are not at the edges anymore:\n')
        outfile.write(str(names))
        outfile.write('\n')
    print(len(names), ' of ',len(X_insteek), ' profiles are dropped because the wet profile is shifted such that A24/A28 are not at the edges anymore. Check ', paths.txt_dropped, ' for the ID codes.')
    
    # Combine Geopandas DataFrames
    gdf_combined = gpd.GeoDataFrame(pd.concat([gdf1, gdf2], ignore_index=True))
    
    # Group based on the column 'profielmet'
    grouped = gdf_combined.groupby('profielmet')
    
    # Iterate through profiles
    for name, group in tqdm(grouped):
        try:
            # Find the point with code 'A22L'
            reference_point = group[group['punttype'] == 'A22L'].geometry.iloc[0]
            afstand_ref     = group[group['punttype'] == 'A22L']['afstand'].values
            afstand_ref     = np.round(afstand_ref,2)
            
            # Estimate the distance (where this is missing) from each point to the point 'A22L'
            afstand = group[group['afstand'].isnull()].geometry.apply(lambda x: x.distance(reference_point))
            afstand = np.round(afstand,2)
            punttype= group[group['afstand'].isnull()]['punttype']
            
            # Left side
                # point with the smallest distance to A22L
                # at that location the distance is subtracted from the distance at A22L
            afstand_rel = afstand.copy()
            afstand_rel.iloc[np.argmin(afstand)] = afstand_ref - afstand.iloc[np.argmin(afstand)]
            punttype.iloc[np.argmin(afstand)]= code_links
            
            # Right side
                # point with the largest distance to A22L
                # at that location the distance is added to the distance at A22L
            afstand_rel.iloc[np.argmax(afstand)] = afstand_ref + afstand.iloc[np.argmax(afstand)]
            punttype.iloc[np.argmax(afstand)]= code_rechts
            
            # Fill in new values ("punttype" and "afstand") in gdf_combined
            gdf_combined.loc[(gdf_combined['profielmet'] == name) & (gdf_combined['afstand'].isnull()), 'punttype'] = punttype
            gdf_combined.loc[(gdf_combined['profielmet'] == name) & (gdf_combined['afstand'].isnull()), 'afstand'] = afstand_rel

        except:
              print('Error profile: ', name)  
              
    # save as shapefile
    gdf_combined = gdf_combined.sort_values(['profielmet','afstand'], ascending=[True,True])
    gdf_combined['profiel'] = gdf_combined[['AHN', 'slibhoogte']].fillna(0).sum(axis=1, numeric_only=True)
    gdf_combined.to_file(paths.shp_points_incl_AHN)
    
def fun_punt_op_afstand(input_data,afstand_nieuw_punt):
    ''' Function to find points at a specific distance from the waterline ("afstand_nieuw_punt") including the corresponding height data''' 
    
    # Input data
    paths           = input_data['paths']
    group           = input_data['group']
    code_insteek    = input_data['code_insteek']
    code_water      = input_data['code_water']
    d_edge          = input_data['d_edge']
    p_insteek       = group[group['punttype'] == code_insteek].geometry
    p_water         = group[group['punttype'] == code_water].geometry
    
    # xy from water edge to profile edge (insteek)
    x = group[(group['punttype'] == code_insteek) | (group['punttype'] == code_water)].geometry.x.values
    y = group[(group['punttype'] == code_insteek) | (group['punttype'] == code_water)].geometry.y.values
    
    # regression line
    fy = np.poly1d(np.polyfit(x, y, 1)) # y = ax + b
    fx = np.poly1d(np.polyfit(y, x, 1)) # x = ay + b
    
    
    
    if abs(fy[1]) < 1.5: # not a vertical line (based on the slope of the regression line)
        
        # Create points from the waterline to the profile edge at an interval of dxy in the x-direction
        # estimate the corresponding y-coordinate based on the function fy
        # estimate the distance from each point to the waterline
        # continue with decreasing interval until the distance deviates max 1 cm from the desired distance (variabel "afstand_nieuw_punt")
        delta_distance = 999
        i_dxy          = 0
        
        # interval step over which to look for points
        c_dir   = np.where(p_insteek.x.values - p_water.x.values>0,1,-1) # direction: >0: left to right, < 0: right to left
        dxy     = [0.5, 0.25, 0.1, 0.05, 0.01, 0.005] * c_dir
        # distance beyond profile edge to look for points
        d_edge  = d_edge * c_dir
        
        while delta_distance > 0.01:
            Xcoords  = np.arange(p_water.x.values,p_insteek.x.values+d_edge,dxy[i_dxy])        
            df       = pd.DataFrame({'X':Xcoords,'Y':fy(Xcoords)})
            gdf      = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.X,df.Y)])                       
            distance = np.round(gdf.geometry.apply(lambda x: x.distance(p_water.iloc[0])),2)
            if (afstand_nieuw_punt > distance.iloc[-1]): 
                ind            = np.argmin(abs(distance-distance.iloc[-1]))
                delta_distance = abs(distance[ind] - distance.iloc[-1])
                distance       = distance[ind]                  
            else:
                ind            = np.argmin(abs(distance-afstand_nieuw_punt))
                delta_distance = abs(distance[ind] - afstand_nieuw_punt)
                distance       = distance[ind]                              
            i_dxy    = i_dxy + 1                               
                
        # find the xy-coordinate of the point that has a distance closest to the desired distance (variabel "afstand_nieuw_punt")
        Xcoord   = Xcoords[ind]
        Ycoord   = fy(Xcoord)        
    else:
        # estimate x-coord from regression (assuming line is vertical)
        delta_distance = 999
        i_dxy          = 0
        
        # interval step over which to look for points
        c_dir   = np.where(p_insteek.y.values - p_water.y.values>0,1,-1) # direction: >0: top to bottom, < 0: bottom to top
        dxy     = [0.5, 0.25, 0.1, 0.05, 0.01, 0.005] * c_dir
        # distance beyond profile edge to look for points
        d_edge  = d_edge * c_dir
        
        while delta_distance > 0.01:
            Ycoords   = np.arange(p_water.y.values,p_insteek.y.values+d_edge,dxy[i_dxy])       
            df       = pd.DataFrame({'X':fx(Ycoords),'Y':Ycoords})
            gdf      = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.X,df.Y)])                       
            distance = np.round(gdf.geometry.apply(lambda x: x.distance(p_water.iloc[0])),2)            
            if (afstand_nieuw_punt > distance.iloc[-1]): 
                ind            = np.argmin(abs(distance-distance.iloc[-1]))
                delta_distance = abs(distance[ind] - distance.iloc[-1])
                distance       = distance[ind]                  
            else:
                ind            = np.argmin(abs(distance-afstand_nieuw_punt))
                delta_distance = abs(distance[ind] - afstand_nieuw_punt)
                distance       = distance[ind]                              
            i_dxy    = i_dxy + 1             
            
        # find the xy-coordinate of the point that has a distance closest to the desired distance (variabel "afstand_nieuw_punt")
        Ycoord   = Ycoords[ind]
        Xcoord   = fx(Ycoord)
    df       = pd.DataFrame({'X':[Xcoord],'Y':[Ycoord],
                             'HO_ID': group['HO_ID'].values[0], 'HO_afstd': group['HO_afstd'].values[0],'HO_cat':group['HO_cat'].values[0]})    
    gdf      = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.X,df.Y)])                
        
    # add AHN
    raster_path    = paths.tif_AHN
    raster         = gdal.Open(raster_path) 
    gt             = raster.GetGeoTransform() 
    raster_band    = raster.GetRasterBand(1) 
    gdf['AHN']     = get_elevation(gdf.geometry.x, gdf.geometry.y, gt, raster_band)
    return gdf

def fun_add_points(gdf, gdf_new, name, a_ref, p_water, oeverkant, punttype):  
    ''' Add new points (gdf_new) to shapefile (gdf)''' 
    
    gdf_new['profielmet'] = name
    gdf_new['profiel']    = gdf_new['AHN']
    gdf_new['punttype']   = punttype
    gdf_new['HO_ID']      = gdf[gdf['profielmet']==name]['HO_ID'].values[0]
    if oeverkant == 'links':
        gdf_new['afstand']    = np.round(a_ref - gdf_new.geometry.iloc[0].distance(p_water),2)
    else:
        gdf_new['afstand']    = np.round(a_ref + gdf_new.geometry.iloc[0].distance(p_water),2)
    gdf                   = gpd.GeoDataFrame(pd.concat([gdf, gdf_new], ignore_index=True))
    gdf                   = gdf.sort_values(['profielmet','afstand'], ascending=[True,True])
    return gdf

def extra_droge_punten(paths):
    ''' Add additional points in the dry part of the profile if the disance from the waterline to the profile edge is large''' 
    
    # Input data
    gdf = gpd.read_file(paths.shp_points_incl_AHN)
    
    # Group points based on the column 'profielmet'
    grouped = gdf.groupby('profielmet')
    
    # Interval in which to add points
    afstanden = [2.5, 5, 7.5] # interval extra point
    d_edge    = 15 # additional distance beyond profile edge (insteek) to be included
    
    # Loop through each profielmet-code
    for name, group in tqdm(grouped):        
            # left side
            p_insteek    = group[group['punttype'] == 'A24'].geometry.iloc[0]
            p_water      = group[group['punttype'] == 'A22L'].geometry.iloc[0]
            a_ref        = group[group['punttype'] == 'A22L']['afstand'].values
            code_extra_dr= 'A31'
            code_rand_dr = 'A23'
            code_insteek = 'A24'
            code_water   = 'A22L'
            afstd        = np.round(p_insteek.distance(p_water),2)
            input_data   = dict({'paths':paths,'d_edge':d_edge,'group':group,'code_insteek':code_insteek,'code_water':code_water})
            # add points between A22L and A24
            for d in range(0,len(afstanden)):              
                if afstd > afstanden[d]:
                    # add point at afstanden[d]
                    afstand_nieuw_punt = afstanden[d]
                    gdf_new            = fun_punt_op_afstand(input_data, afstand_nieuw_punt)
                    gdf                = fun_add_points(gdf, gdf_new, name, a_ref, p_water, 'links', code_extra_dr) # toevoegen aan shapefile 
                
                
            # add points beyond A24
            AHN_prev = group[group['punttype']==code_insteek]['profiel'].values
            for d in range(0,len(afstanden)):  
                if (afstanden[d] > afstd) & (afstanden[d] < afstd + d_edge):
                    afstand_nieuw_punt = afstanden[d]
                    gdf_new            = fun_punt_op_afstand(input_data, afstand_nieuw_punt)
                    if gdf_new['AHN'].values > AHN_prev:                    
                        gdf                = fun_add_points(gdf, gdf_new, name, a_ref, p_water, 'links', code_rand_dr) # toevoegen aan shapefile                         
                        AHN_prev           = gdf_new['AHN'].values
                    else:
                        break
            
            # right side
            p_insteek    = group[group['punttype'] == 'A28'].geometry.iloc[0]
            p_water      = group[group['punttype'] == 'A22R'].geometry.iloc[0]
            a_ref        = group[group['punttype'] == 'A22R']['afstand'].values
            code_extra_dr= 'A32'
            code_rand_dr = 'A29'
            code_insteek = 'A28'
            code_water   = 'A22R'        
            afstd        = np.round(p_insteek.distance(p_water),2)        
            input_data   = dict({'paths':paths,'d_edge':d_edge,'group':group,'code_insteek':code_insteek,'code_water':code_water})
            
            # add points between A22R and A28
            for d in range(0,len(afstanden)):                
                if afstd > afstanden[d]:
                    # add point at afstanden[d]
                    afstand_nieuw_punt = afstanden[d]
                    gdf_new            = fun_punt_op_afstand(input_data, afstand_nieuw_punt)
                    gdf                = fun_add_points(gdf, gdf_new, name, a_ref, p_water, 'rechts',code_extra_dr ) # toevoegen aan shapefile 
            
            # add points beyond A28
            AHN_prev = group[group['punttype']==code_insteek]['profiel'].values
            for d in range(0,len(afstanden)):                        
                if (afstanden[d] > afstd) & (afstanden[d] < afstd + d_edge):
                    afstand_nieuw_punt = afstanden[d]
                    gdf_new            = fun_punt_op_afstand(input_data, afstand_nieuw_punt)                    
                    if gdf_new['AHN'].values > AHN_prev:                    
                        gdf                = fun_add_points(gdf, gdf_new, name, a_ref, p_water, 'rechts', code_rand_dr) # toevoegen aan shapefile                                         
                        AHN_prev           = gdf_new['AHN'].values
                    else:
                        break
        
        
    # save file
    gdf.to_file(paths.shp_points_AHN_all)    
    

        