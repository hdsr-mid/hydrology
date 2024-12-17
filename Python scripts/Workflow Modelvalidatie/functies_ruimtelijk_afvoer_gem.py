# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:20:57 2024

@author: PetraH

gebaseerd op scripts ontvangen voor het berekenen van debietsstatistieken
"""
import os
import numpy as np
import geopandas as gpd
import hkvsobekpy as hkv
import pandas as pd
import bokeh.plotting 
from bokeh.models import ColumnDataSource, GeoJSONDataSource, LinearColorMapper, Range1d, CustomJS
from bokeh.models import HoverTool, Whisker, VBar, TapTool, ColorBar
from bokeh.layouts import layout
from rasterstats import zonal_stats
from osgeo import gdal

import warnings
warnings.filterwarnings("ignore")

bokeh.plotting.output_notebook()


vmin   = -2
vmax   = 2

# voor de boxplot en histogram
dmin  = -1000
dmax  = 5000
nbins = dmax-dmin

def bokeh_boxplot(shp_comp,type_OW):
    if type_OW=='oost':
        col_Q_WB = 'Q_WB_O'        
    else:
        col_Q_WB = 'Q_WB_W'        
    
    # create dataframe
    df1 = pd.DataFrame({'kind': ['Water balans']*len(shp_comp),'value': shp_comp[col_Q_WB]})
    df2 = pd.DataFrame({'kind': ['HYDROMEDAH']*len(shp_comp),'value': shp_comp['Q_HYDRMDH']})
    df  = pd.concat([df1, df2],axis=0, join="inner",ignore_index=False)     
    kinds = df.kind.unique()
    
    # compute quantiles
    qs = df.groupby("kind").value.quantile([0.25, 0.5, 0.75])
    qs = qs.unstack().reset_index()
    qs.columns = ["kind", "q1", "q2", "q3"]
    df = pd.merge(df, qs, on="kind", how="left")

    # compute IQR outlier bounds
    iqr = df.q3 - df.q1
    df["upper"] = df.q3 + 1.5*iqr
    df["lower"] = df.q1 - 1.5*iqr

    source = ColumnDataSource(df)

    p = bokeh.plotting.figure(x_range=kinds, tools="", toolbar_location=None,y_axis_label="Q [mm/yr]", height=350, width=350)

    # outlier range
    whisker = Whisker(base="kind", upper="upper", lower="lower", source=source)
    whisker.upper_head.size = whisker.lower_head.size = 20
    p.add_layout(whisker)

    # quantile boxes
    p.vbar("kind", 0.7, "q2", "q3", source=source, line_color="black")
    p.vbar("kind", 0.7, "q1", "q2", source=source, line_color="black")

    # outliers
    outliers = df[~df.value.between(df.lower, df.upper)]
    p.scatter("kind", "value", source=outliers, size=6, color="black", alpha=1)

    p.xgrid.grid_line_color = None
    p.axis.major_label_text_font_size="14px"
    p.axis.axis_label_text_font_size="12px"
    p.y_range = Range1d(dmin,dmax) 

    return p
    
def bokeh_fig(paths, shp_afvoer, shp_comp, type_OW):  
    if type_OW=='oost':
        col_kwel = 'Kwel_oost'
        col_Q_WB = 'Q_WB_O'
        metric   = 'dQ_O'        
        fightml  = paths.fightml.replace('.html','_oost.html')
    else:
        col_kwel = 'Kwel_west'
        col_Q_WB = 'Q_WB_W'
        metric   = 'dQ_W'        
        fightml  = paths.fightml.replace('.html','_west.html')
    
    # remove polygons out of region
    shp_comp_sel = shp_comp.copy()
    for col in ['Psat','Esat','Q_HYDRMDH']:
        shp_comp_sel[np.isnan(shp_comp[col_Q_WB])][col] = np.nan 
    
    # -----------------------------------------------
    # Get data
    shp_comp_sel['FID'] = shp_comp_sel.reset_index().index.values.astype(str)
    geo_source = GeoJSONDataSource(geojson=shp_afvoer.to_json())
    map_source = GeoJSONDataSource(geojson=shp_comp_sel.to_json())
    
    # -----------------------------------------------
    # --------- Plot map --------- 
    p1 = bokeh.plotting.figure(title='Afvoer (laterals) gemiddeld: model (HYDROMEDAH) vs. referentie (satelliet)', height=350, width=820)
    color = LinearColorMapper(palette = 'Turbo256', low = vmin, high = vmax)
    map_poly  = p1.patches(fill_alpha=1,
              fill_color={'field': metric, 'transform': color},
              line_color='black', line_width=0.5, source=map_source)
    color_bar = ColorBar(color_mapper=color,title='(model - referentie)/referentie [-]')
    p1.add_layout(color_bar, 'right')
    tooltips = [('ID', '@CODE'),
                ('Label', '@'+metric+'{0.}'),
                ('HYDROMEDAH', '@Q_HYDRMDH{0.}'+'mm/yr'),
                ('Water balans', '@'+col_Q_WB+'{0.}'+'mm/yr')]
    hover = HoverTool(renderers=[map_poly], tooltips=tooltips) 
    p1.add_tools(hover)  
    
    # --------- boxplot --------- 
    p2 = bokeh_boxplot(shp_comp_sel,type_OW)
    
    # --------- bar plot --------- 
    cols     = ['Psat',col_kwel,'Esat',col_Q_WB,'Q_HYDRMDH']
    shp_comp_sel['Esat']         = -shp_comp_sel['Esat']
    shp_comp_sel[col_Q_WB]       = -shp_comp_sel[col_Q_WB]
    shp_comp_sel['Q_HYDRMDH']    = -shp_comp_sel['Q_HYDRMDH']
    src_dict = fun_create_dict(shp_comp_sel,cols)
    ID_0     = shp_comp_sel.iloc[0]['FID']
    src      = ColumnDataSource(data=src_dict[ID_0])        
    
    p3 = bokeh.plotting.figure(x_range=cols,title='Waterbalans ' + shp_comp_sel.iloc[0]['CODE'], height=350, width=820,y_axis_label="Flux [mm/yr]")
    glyph = VBar(x='index', top='data', bottom=0, width=0.5, fill_color='powderblue')
    p3.add_glyph(src, glyph)
    
    #create the TapTool, add it to the bar figure, and enable it on initialization
    tap = TapTool(renderers=[map_poly])
    p1.add_tools(tap)
    p1.toolbar.active_tap = tap
    
    cb = CustomJS(args=dict(map_source=map_source,src_dict=src_dict,src=src,p3=p3,glyph=glyph)
                  ,code='''
                  var sel_bar_i = map_source.selected.indices[0]
                  var line_id = map_source.data['FID'][sel_bar_i]
                  src.data = src_dict[line_id]
                  p3.title.text= 'Waterbalans ' + map_source.data['CODE'][sel_bar_i]
                  src.change.emit()
                  p3.change.emit()
                  ''')
    
    
    #apply this callback to trigger whenever the indices of the selected glyph change
    map_source.selected.js_on_change('indices',cb)             
    
    grid = layout([[p1,p2],[p3,]])
    bokeh.plotting.show(grid, notebook_handle=True)
    bokeh.plotting.save(grid,fightml)    
    
def add_statistic(output, data, threshold=None, percentile=None, name=None):
    # functie om statistieken toe te voegen aan de shapefile
    data = abs(data)
    
    statcol = data.mean(axis=0, skipna=True)
    statcol.name = name
    output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)     
    
    return(output)

def data_2_shp(shp, df):
    ID      = df['ID'].values    
    nmax    = len(shp)
    shp_dat = np.zeros((nmax, 1))-999
    for i in range(0,len(df)):
        ind = np.where(ID[i]==shp['ID'].values)
        if len(ind)==1: 
            shp_dat[ind[0]]= df['Q_Sobek'].values[i]         
        else:
            print(ID[i], "uitzondering 2, weggelaten...")
    shp['Q_Sobek']  = shp_dat
    return shp

def lat_data(paths):
    data_file = open(paths.dat_lat,'r')
    ID   = []
    area = []
    rain = []
    sep  = []
    for line in data_file:
        data = line.split()
        ID.append(str(data[2].replace("'","")))
        area.append(float(data[-2])) #m2
        rain.append(float(data[11])) #mm/s
        sep.append(float(data[15])) #mm/s
    
    rain = np.array(rain)
    sep  = np.array(sep)
    Qlat = (rain + sep)/ 1000 * area # m3/s
    df_Qlat = pd.DataFrame({'ID': ID,'Q_Sobek': np.round(Qlat,2)})
    return df_Qlat

def fun_create_dict(shp,cols):
    # Cerate timeseries dictionary
    src_dict          = {}
    
    # Create dictionary
    for point in range(0,len(shp)):    
        FID           = shp.iloc[point]['FID']        
        data          = shp.iloc[point][cols].T.reset_index()
        data.columns  = ['index','data']
        src_dict[FID] = dict(data)
    
    return src_dict

    
def get_raster_pixel_size(path):
    # load raster
    raster = gdal.Open(path)
    
    # get pizel size
    gt =raster.GetGeoTransform()
    
    # close raster
    raster = None
    
    return abs(gt[1])

def get_kwel(paths, shp_afvoer):
    # copy to avoid overwriting
    shp = shp_afvoer.copy()
    
    # pixel size
    dx1 = get_raster_pixel_size(paths.raster_kwelO)
    dx2 = get_raster_pixel_size(paths.raster_kwelW)
    
    # Bereken median kwel
    statsO          = zonal_stats(shp,paths.raster_kwelO,stats=['median','count'])
    statsW          = zonal_stats(shp,paths.raster_kwelW,stats=['median','count'])
    shp             = pd.concat((shp, pd.DataFrame(statsO)), axis=1).rename(columns ={'median':'Kwel_oost','count':'Aoost'})
    shp             = pd.concat((shp, pd.DataFrame(statsW)), axis=1).rename(columns ={'median':'Kwel_west','count':'Awest'})
    
    # conversions
    shp['Area']      = shp.area
    shp['Aoost']     = shp['Aoost'] * dx1*dx1
    shp['Awest']     = shp['Awest'] * dx2*dx2
    shp['ratio_oost']= shp['Aoost'] / shp.area
    shp['ratio_west']= shp['Awest'] / shp.area
    shp['Kwel_oost'] = np.where(shp['ratio_oost'] > 0.75, shp['Kwel_oost']*365, -99999)
    shp['Kwel_west'] = np.where(shp['ratio_west'] > 0.75, shp['Kwel_west']*365, -99999)
    shp['Kwel_oost'] = np.where((shp['ratio_oost'] <= 0.75) & (shp['ratio_west'] <= 0.75), -99999, shp['Kwel_oost'])    
    shp['Kwel_west'] = np.where((shp['ratio_oost'] <= 0.75) & (shp['ratio_west'] <= 0.75), -99999, shp['Kwel_west'])    
    shp              = shp.drop(['Area','Aoost', 'Awest','ratio_oost','ratio_west'], axis=1)
    
    # save temp file
    shp.to_file(paths.shp_kwel) 
    
    return shp.replace(-99999,np.nan)

def main(paths, months, years):
    ''' Satellite data: FEWS-WIS'''
    # input data
    shp_afvoer = gpd.read_file(paths.shp_afvoer)
    shp_afwatr = gpd.read_file(paths.shp_afwatr)
    shp_peil   = gpd.read_file(paths.shp_peil).rename(columns ={'GPGIDENT':'CODE'})
    df_P       = pd.read_csv(paths.df_P, delimiter=',')
    df_E       = pd.read_csv(paths.df_E, delimiter=',')
    
    # format data
    shp_afvoer.crs ='EPSG:28992'
    shp_afwatr.crs ='EPSG:28992'    
    shp_peil.crs   ='EPSG:28992'    
    df_P   = df_P.drop(0).drop(columns='GMT+1').astype(float)
    df_E   = df_E.drop(0).drop(columns='GMT+1').astype(float)
    
    # calculate mean
    df_P  = df_P.replace(-9999,np.nan).mean(axis=0, skipna=True).reset_index()
    df_E  = df_E.replace(-9999,np.nan).mean(axis=0, skipna=True).reset_index()    
    df_P.columns = ['CODE','Psat']
    df_E.columns = ['CODE','Esat']
    
    # add values to shapefile
    shp_sat_peil   = shp_peil[['CODE','geometry']].merge(df_P[['CODE','Psat']], on='CODE', how='left')
    shp_sat_peil   = shp_sat_peil.merge(df_E[['CODE','Esat']], on='CODE', how='left')
    
    # weighted mean op afvoergebiedsniveau
    cols                     = ['Psat','Esat']
    shp_sat_afvoer           = shp_afvoer.copy()
    shp_sat_afvoer['A_afvr'] = shp_sat_afvoer.area    
    shp_sat_peil['C_peil']   = shp_sat_peil['CODE']
    shp_sat_peil             = shp_sat_peil[['C_peil','Psat','Esat','geometry']]
    df_is                    = gpd.overlay(shp_sat_peil, shp_sat_afvoer[['CODE','A_afvr','geometry']], how='intersection') 
    df_is['A_peil']          = df_is.area
    for c in cols: 
        df_is[c]       = df_is[c] * df_is['A_peil']/df_is['A_afvr']    
        df_c           = df_is.groupby(['CODE'])[[c]].sum() .reset_index()
        shp_sat_afvoer = shp_sat_afvoer.merge(df_c[['CODE',c]], on='CODE', how='left')
    
    ''' Satellite data: Kwel'''
    shp_sat_afvoer = get_kwel(paths, shp_sat_afvoer)
    
    # conversions
    shp_sat_afvoer['Kwel_oost'] = np.round(shp_sat_afvoer['Kwel_oost'])
    shp_sat_afvoer['Kwel_west'] = np.round(shp_sat_afvoer['Kwel_west'])
    shp_sat_afvoer['Psat']      = np.round(shp_sat_afvoer['Psat'])
    shp_sat_afvoer['Esat']      = np.round(shp_sat_afvoer['Esat'])
    shp_sat_afvoer['Q_WB_O']    = np.round((shp_sat_afvoer['Psat'] + shp_sat_afvoer['Kwel_oost'] - shp_sat_afvoer['Esat']))
    shp_sat_afvoer['Q_WB_W']    = np.round((shp_sat_afvoer['Psat'] + shp_sat_afvoer['Kwel_west'] - shp_sat_afvoer['Esat']))
    
    # shapefile opslaan
    shp_sat_afvoer.to_file(paths.shp_sat)   
    
    '''HYDROMEDAH'''
    # shapefiles inlezen
    shp_line = gpd.read_file(paths.shp_RchSeg)
    shp_line.crs ='EPSG:28992'
    shp_line.columns = [col.strip() for col in shp_line.columns]
    shp_line.index = shp_line['ID']
    shp_line.drop(['NAME', "TYPE", 'PARENTID', 'ID_FROM', 'ID_TO', 'USERID'], axis=1, inplace=True)
    shp_node  = gpd.read_file(paths.shp_node)
    shp_node.crs ='EPSG:28992'
    shp_node.columns = [col.strip() for col in shp_node.columns]
    
    # Data inlezen 
    # (verwerkt in Sobek HIS bestanden)    
    cases = [1,2,3,4,5,6,7,8,9,10,11,12]#
    subhislist = []
    for case in cases:
        pad_to_his = os.path.join(paths.dat_Sobek, str(case), 'reachseg.his')
        reachseg = hkv.read_his.ReadMetadata(pad_to_his)    
        reach_segments = reachseg.DataFrame()
        reach_segments = reach_segments.iloc[32:,:]
        locs= reachseg.GetLocations()
        params= reachseg.GetParameters()
        reach_segments = reach_segments.loc[:,params[params.index('Discharge mean(m³/s)')]]
        reach_segments.head()
        subhislist.append(reach_segments)
    all_data = pd.concat(subhislist)
    all_data['MONTH'] = all_data.index.month
    all_data['YEAR'] = all_data.index.year
    
    # Statistiek berekenen
    sub_data = all_data[all_data['MONTH'].isin(months)]
    sub_data = all_data[all_data['YEAR'].isin(years)]
    shp_HYDROMEDAH = add_statistic(shp_line, sub_data, name='Q_HYDRMDH')
    
    # Data overzetten naar punten (spatial join)
    shp_HYDROMEDAH = shp_node.sjoin_nearest(shp_HYDROMEDAH, how="left", distance_col="Distances")
    
    # check distance
    ind = np.where(shp_HYDROMEDAH['Distances'].values>0.01)[0]
    if len(ind)>0: 
        print('Check file for nodes with distance > 0.01 m! : ')
        print(paths.shp_HYDROMEDAH)
        shp_HYDROMEDAH.to_file(paths.shp_HYDROMEDAH)
    
    # convert based on area afwateringsgebied (m3/s -> mm/yr)
    shp_afwatr['area']              = shp_afwatr.area
    shp_HYDROMEDAH['CODE']          = [i.replace('drain','') for i in shp_HYDROMEDAH.ID_left.values]    
    shp_HYDROMEDAH                  = shp_HYDROMEDAH.merge(shp_afwatr[['CODE','area']], on='CODE', how='left')
    shp_HYDROMEDAH['Q_HYDRMDH']     = np.round(shp_HYDROMEDAH['Q_HYDRMDH']/shp_HYDROMEDAH['area'] * 1e3 * 3600*24*365,1)
    
    # shapefile opslaan
    shp_HYDROMEDAH.to_file(paths.shp_HYDROMEDAH)    
   
    # Data overzetten naar shp_afvoer (spatial join)    
    shp_HYDROMEDAH     = shp_HYDROMEDAH.rename(columns ={'ID_left':'ID1','ID_right':'ID2'})[['ID1','ID2','Q_HYDRMDH','geometry']]
    shp_HYDROMEDAH     = gpd.sjoin(shp_HYDROMEDAH, shp_afvoer[['CODE','geometry']], how="left")
    shp_HYDROMEDAH     = shp_HYDROMEDAH.groupby('CODE')['Q_HYDRMDH'].agg(['mean']).reset_index()
    shp_HYDROMEDAH.columns = ['CODE','Q_HYDRMDH']
    
    '''Samenvoegen & verschil berekenen'''
    shp_comp = shp_sat_afvoer.merge(shp_HYDROMEDAH, on='CODE', how='left')
    shp_comp['dQ_O'] = np.round((shp_comp['Q_HYDRMDH'] - shp_comp['Q_WB_O'])/shp_comp['Q_WB_O'])
    shp_comp['dQ_W'] = np.round((shp_comp['Q_HYDRMDH'] - shp_comp['Q_WB_W'])/shp_comp['Q_WB_W'])
     
    # save file    
    shp_comp.to_file(paths.shp_comp)    
    # shp_comp = gpd.read_file(paths.shp_comp)
    # Data plotten
    bokeh_fig(paths, shp_afvoer,shp_comp, 'oost')
    bokeh_fig(paths, shp_afvoer,shp_comp, 'west')
    
    