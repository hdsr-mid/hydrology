"""
Generate DFM_lateral_sources.bc from SIMGRO output (sw_dtsw or sw_dtgw)
Base script: "1_1_generate_lateralsources_bc_from_simgro_output.py"
Location: reken05 - "D:\Bovenregionale_stresstest_ARKNZK\0_scripts"

Author:         Petra Hulsman (HDSR)
Last update:    13/04/2026
environment: env_lats_Petra
"""


import pandas as pd
import numpy as np
import geopandas as gpd
from functions import read_bda

import bokeh.plotting 
from bokeh.layouts import layout
from bokeh.models import ColumnDataSource, GeoJSONDataSource, LinearColorMapper, ColorBar
from bokeh.models import HoverTool, TapTool, CustomJS, Line

class paths():
    shp_afvoergebieden      = 'D:/Petra/lats_fig/Input/Afvoergebieden_30032026.gpkg'
    shp_afwateringsgebieden = 'D:/Petra/lats_fig/Input/Afwateringseenheden_20240606.gpkg'
    simgro_path             = r'D:\hydromedah_P23006_KALIB_V4B_LONGPERIOD_20241024_20110101_20240101_250m\output\metaswap'
    outdir                  = 'D:/Petra/lats_fig/Output/'

def fun_create_dict(df):
    # Create timeseries dictionary
    src_dict          = {}
    var               = ['x','y']
    
    for c in df.columns:    
         
        data  = [df.index.values, df[c].values]
                     
        # Create dictionary
        src_dict['AW' + c.zfill(4)] = dict(zip(var, data))
        
    
    return src_dict
    
if __name__ == "__main__":

    print('Start')
    # Input data
    shp_afvoergebied = gpd.read_file(paths.shp_afvoergebieden)
    shp_afwateringsgebieden = gpd.read_file(paths.shp_afwateringsgebieden)
    
    # for ig in range(0,10):
    for ig in range(0,len(shp_afvoergebied)):
        print(ig,len(shp_afvoergebied))
        # select region
        gebied = shp_afvoergebied.iloc[ig]
        CODE   = gebied.CODE
        NAAM   = gebied.NAAM
        shp_afw_sel = shp_afwateringsgebieden.clip(gebied.geometry)
        shp_afv_sel = shp_afvoergebied.drop([ig])
        sel_loc = [str(int(i.replace('AW',''))) for i in shp_afw_sel.CODE.values]
        
        # load laterals from SIMGRO output        
        df = read_bda(simgro_path=paths.simgro_path, 
                      sel_loc = sel_loc,
            filename='sw_dtsw')
        
        # T1
        months = [pd.to_datetime(t).month for t in df.index.values]
        summer = np.array([i for i in range(0,len(months)) if months[i] in np.array([5,6,7,8,9])])
        df_sum = df.iloc[summer]
        year   = [pd.to_datetime(t).year for t in df_sum.index.values]
        df_T1  = df_sum.groupby(year).max().median() # m3/s
        shp_afw_sel['T1'] = -999.0
        for idx in shp_afw_sel.index.values:
            g = shp_afw_sel.loc[idx,'CODE']
            try:
                shp_afw_sel.loc[idx,'T1'] = float(df_T1[str(int(g.replace('AW','')))])
            except:
                shp_afw_sel.loc[idx,'T1'] = np.nan
        
        # dataframe to dict
        src_dict = fun_create_dict(df)
        tmax     = len(df.index.values)
        for aw in shp_afw_sel.CODE.values: 
            if aw in list(src_dict.keys()):
                src_dict[aw]['T1'] = shp_afw_sel.loc[shp_afw_sel['CODE']==aw]['T1'].values*np.ones(tmax,)
            
        # plot - basemap
        geo_source1  = GeoJSONDataSource(geojson=shp_afv_sel.to_json())
        geo_source2  = GeoJSONDataSource(geojson=shp_afw_sel.to_json())
        map_source  = ColumnDataSource(shp_afw_sel)
        
        p1 = bokeh.plotting.figure(title='Afvoer (laterals)', height=350, width=820, match_aspect=True)
        color = LinearColorMapper(palette = 'Turbo256', low = -2, high = 2)
        # map_poly1 = p1.patches(fill_alpha=1, fill_color = 'white',line_color='black', line_width=0.5, source=geo_source1)
        map_poly2 = p1.patches(fill_alpha=1, fill_color={'field': 'T1', 'transform': color},line_color='black', line_width=0.5, source=geo_source2)
        color_bar = ColorBar(color_mapper=color,title='Lateral flow [m3/s]')
        p1.add_layout(color_bar, 'right')
        # p1.add_tools(HoverTool(renderers=[map_poly1], tooltips=[('ID', '@CODE')]) )  
        p1.add_tools(HoverTool(renderers=[map_poly2], tooltips=[('ID', '@CODE'),('T1', '@T1{0.01}' + ' m3/s')]) )  
        
        # plot - timeseries
        # Tapping...
        #start with initial state for columndatasource
        ID_0  = list(src_dict.keys())[0]
        src   = ColumnDataSource(data=src_dict[ID_0])
        p3    = bokeh.plotting.figure(title=ID_0, height=300, width=1400, x_axis_type='datetime', x_axis_label="time",y_axis_label="Q [m3/s]") # , y_range=(0,1)
        glyph = Line(x="x", y='y', line_color="black")
        p3.add_glyph(src, glyph)
        glyph = Line(x="x", y='T1', line_color="red")
        p3.add_glyph(src, glyph)
        
        #create the TapTool, add it to the bar figure, and enable it on initialization
        tap = TapTool(renderers=[map_poly2])
        p1.add_tools(tap)
        p1.toolbar.active_tap = tap
        
        cb = CustomJS(args=dict(map_source=geo_source2,src_dict=src_dict,src=src,p3=p3,glyph=glyph)
                      ,code='''
                      var sel_bar_i = map_source.selected.indices[0]
                      var area_id = map_source.data['CODE'][sel_bar_i]
                      p3.title.text= area_id
                      p3.change.emit()
                      src.data = src_dict[area_id]
                      src.change.emit()
                      ''')
        
        # #apply this callback to trigger whenever the indices of the selected glyph change
        geo_source2.selected.js_on_change('indices',cb)             
        
        grid = layout([[p1],[p3,],])
        # bokeh.plotting.show(grid, notebook_handle=True)
        bokeh.plotting.save(grid,paths.outdir + str(ig) + 'fig_Qlats_' + CODE + '_' + NAAM + '.html')
        
            
        
    print('Done')
