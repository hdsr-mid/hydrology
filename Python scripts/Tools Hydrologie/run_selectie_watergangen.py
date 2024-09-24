# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 12:15:19 2024

@author: phulsman
"""

import os
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
from typing import List
from tqdm import tqdm

class paths():
    root = os.getcwd()
    
    shp_afvoergbd   = os.path.join(root,'GIS_input','br_vs_afvoergebieden.shp')
    shp_watergangen = os.path.join(root,'GIS_input','br_vs_hydroobject.shp')
    shp_RUPROF      = os.path.join(root,'GIS_input','Basis_Hydraulische_Toetsing_Hydro_PG_SOBEK_with_Dis_ZOMER_20240715.shp')
    shp_boezem      = os.path.join(root,'GIS_input','br_vs_boezemwatersysteem.shp')
    
    shp_out         = os.path.join(root,'GIS_output','resultaten.gpkg')
    
def group_lines(shp):
    if 'group' not in shp.columns:
        shp['group'] = 0
    starts = [str(row.geometry.coords[0]) for index, row in shp.iterrows()]
    ends   = [str(row.geometry.coords[-1]) for index, row in shp.iterrows()]
    
    before = shp.group.values.copy()
    for i in tqdm(range(0,len(shp))):
        cond = ((starts[i] == np.array(starts)) | (starts[i] == np.array(ends)) | (ends[i] == np.array(starts)) | (ends[i] == np.array(ends)))
        idx  = np.where(cond)[0]
        
        group_values = np.unique(shp.iloc[idx]['group'].values)
        group_values = group_values[group_values!=0]
        if len(group_values)==0:
            shp.loc[cond, 'group']   = i
        else:
            shp.loc[cond, 'group'] = int(np.min(group_values))
    
    after  = shp.group.values.copy()
    changes= len(np.where(before != after)[0])
    
    # i = int(np.where(shp.CODE == 'H071905')[0])
    # cond = ((starts[i] == np.array(starts)) | (starts[i] == np.array(ends)) | (ends[i] == np.array(starts)) | (ends[i] == np.array(ends)))
    # idx  = np.where(cond)[0]
    # print(idx)
    # print(shp.iloc[idx])
    
    return shp, changes

if __name__ == "__main__":
        # thresholds
        width_percentile = 0.5
        Q_threshold      = 0.5
        length_perc      = 0.5
        
        # input data
        col_watergangen = ['CODE','CATEGORIEO','BREEDTE','WS_MODEL','geometry']
        col_afvoergbd   = ['CODE','geometry']
        col_RUPROF      = ['hydro_br','Value','geometry']
        col_boezem      = ['CODE']
        shp_watergangen = gpd.read_file(paths.shp_watergangen)[col_watergangen]  
        shp_afvoergbd   = gpd.read_file(paths.shp_afvoergbd)[col_afvoergbd]
        shp_RUPROF      = gpd.read_file(paths.shp_RUPROF)[col_RUPROF].rename(columns ={'Value':'Q'})
        shp_boezem      = gpd.read_file(paths.shp_boezem)[col_boezem]
        shp_watergangen.to_file(paths.shp_out, driver="GPKG",layer = 'orig')
        
        # add width information
        shp_watergangen['centroid'] = shp_watergangen.geometry.centroid #Create a centroid point column
        shp_watergangen['polygeom'] = shp_watergangen.geometry #Save the polygon geometry to switch back to after the join
        shp_watergangen             = shp_watergangen.set_geometry('centroid')
        shp_watergangen             = shp_watergangen.sjoin(shp_afvoergbd[['CODE','geometry']].rename(columns ={'CODE':'CODE_afv'}), how="left")
        shp_watergangen['geometry'] = shp_watergangen['polygeom']
        shp_watergangen             = shp_watergangen.set_geometry('geometry').drop(columns=['polygeom', 'centroid'])
        shp_watergangen_sel = shp_watergangen[shp_watergangen['WS_MODEL'] == 'j']        
        p = shp_watergangen_sel[['CODE_afv','BREEDTE']].groupby('CODE_afv').quantile(width_percentile)[['BREEDTE']]
        shp_watergangen     = shp_watergangen.set_index('CODE_afv').join(p, rsuffix='_stat').reset_index()
        
        # add Q data
        shp_RUPROF      = shp_RUPROF[(shp_RUPROF.hydro_br.values.astype(str)!='nan')].rename(columns ={'hydro_br':'CODE'})
        shp_RUPROF      = shp_RUPROF.dissolve(by='CODE', aggfunc='mean').reset_index().drop('geometry', axis=1)        
        shp_watergangen = shp_watergangen.set_index('CODE').join(shp_RUPROF.set_index('CODE'), rsuffix='_RUPROF').reset_index()
        
        # add boezem info
        shp_boezem['BOEZEM'] = 1
        shp_watergangen      = shp_watergangen.set_index('CODE').join(shp_boezem.set_index('CODE')).reset_index()
        
        # conditions
        cond0 = (shp_watergangen['WS_MODEL'] == 'j')
        cond1 = (shp_watergangen['BREEDTE']> shp_watergangen['BREEDTE_stat'])
        cond2 = (shp_watergangen['CATEGORIEO'] == 98)
        cond3 = (shp_watergangen['Q'] > Q_threshold)
        cond4 = (shp_watergangen['BOEZEM'] == 1)
        shp_watergangen['sel'] = np.where(cond2 | (cond0 & cond1) | (cond0 & cond3) | cond4 ,1,0)
        shp_watergangen.to_file(paths.shp_out, driver="GPKG",layer = 'incl selection criteria')
        
        # select waterways
        shp_watergangen_sel = shp_watergangen[shp_watergangen['sel']==1].drop(columns='index_right')
        
        # add group number
        changes = len(shp_watergangen_sel)
        while changes != 0:
            print(changes)
            shp_watergangen_sel, changes = group_lines(shp_watergangen_sel)            
        shp_watergangen_sel.to_file(paths.shp_out, driver="GPKG",layer = 'incl group nr')
        
        # add length info
        shp_watergangen_sel['lengte'] = shp_watergangen_sel.length
        p                             = shp_watergangen_sel.drop(columns=['geometry']).groupby('group')[['lengte']].sum()
        shp_watergangen_sel           = shp_watergangen_sel.drop(columns=['lengte']).set_index('group').join(p).reset_index()
        
        # select waterways
        shp_watergangen_cleaned = shp_watergangen_sel[shp_watergangen_sel['lengte'] > np.quantile(shp_watergangen_sel['lengte'], length_perc)]
        shp_watergangen_cleaned.to_file(paths.shp_out, driver="GPKG",layer = 'cleaned')
        
        # dissolve
        shp_watergangen_sel = shp_watergangen_sel.dissolve(by='group').reset_index() # HO codes kloppen niet meer door deze stap!!!
        shp_watergangen_sel[['geometry']].to_file(paths.shp_out, driver="GPKG",layer = 'dissolved')
        







    # from shapely import MultiLineString, line_merge
    # llout = line_merge(shp_watergangen_sel.geometry)
    # shp_watergangen_sel2 = group_lines(shp_watergangen_sel.explode())
    # shp_watergangen_sel2 = shp_watergangen_sel2.dissolve(by='group')
    # print(len(shp_watergangen_sel2))
    
    
    # from shapely import geometry
    # from shapely.ops import linemerge
    # multi_line = geometry.MultiLineString(list(shp_watergangen_sel.geometry))
    # merged_line = linemerge(multi_line)
    
    # from shapely.geometry import MultiLineString, Point
    # s = gpd.GeoSeries([
    # MultiLineString([[(0, 2), (0, 10)], [(0, 10), (5, 10)]]),
    # MultiLineString([[(0, 2), (0, 10)], [(0, 11), (5, 10)]]),        
    # MultiLineString([[(0, 0), (1, 0)], [(0, 0), (3, 0)]]),])
    # ss = linemerge(s)
    # un = unary_union(shp_watergangen_sel)
    # from scipy.sparse.csgraph import connected_components
    # overlap_matrix = shp_watergangen_sel.apply(lambda x: shp_watergangen_sel.overlaps(x)).values.astype(int)
    # n, ids = connected_components(shp_watergangen_sel)
    