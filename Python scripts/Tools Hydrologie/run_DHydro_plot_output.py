# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:57:52 2024

@author: PetraH

Plot verschillende onderdelen van DHydro:
    - .pliz         fixed weir (hoogtelijnen)
    - .bc           boundary condition
    - _his.nc       output time series 
    - _map.nc       output gridded
    
virtual environment: dsd_env

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import contextily as ctx
import xugrid as xu

if __name__ == "__main__":
    
    
    
    ''' .pliz file (fixed weir/hoogtelijnen) ''' 
    fn = r'DATA/objects.pli_fxw.pliz'
    polyfile_object = hcdfm.PolyFile(fn)
    gdf_polyfile = dfmt.PolyFile_to_geodataframe_linestrings(polyfile_object,crs='EPSG:28992')    
    ax = gdf_polyfile.plot()
    ax.set_aspect('equal')
    plt.savefig(r'figs_output/fig_pliz.png')
    plt.close()
    
    ''' .bc file (boundary condition) ''' 
    fn  = r'DATA/rhine.bc'
    var = 'lateral_discharge'
    forcingmodel_object = hcdfm.ForcingModel(fn)
    fig, ax = plt.subplots()
    for iFO, forcingobj in enumerate(forcingmodel_object.forcing):
        forcing_xr = dfmt.forcinglike_to_Dataset(forcingobj, convertnan=True)
        forcing_xr[var].plot(ax=ax, label=forcing_xr[var].attrs['locationname'], linewidth=0.8)
    ax.legend(loc=1)    
    plt.savefig(r'figs_output/fig_bc.png')
    plt.close()
    
    ''' _his.nc file '''
    fn        = r'DATA/DFM_his.nc'
    var_dat   = ['weirgen_discharge','lateral_realized_discharge_average','dambreak_discharge'] 
    ds_his = xr.open_mfdataset(fn, preprocess=dfmt.preprocess_hisnc)
    # print(ds_his)
    # print(list(ds_his.variables.keys()))
    # print(list(ds_his._coord_names))
    for i in range(0,len(var_dat)):
        ds_his_sel = ds_his[var_dat[i]]
        var_loc    = list(ds_his_sel.coords)[1]
        nmax       = len(ds_his_sel[var_loc])
        if nmax>1:
            ds_his_sel = eval('ds_his_sel.isel('+var_loc+'=np.arange(0,10))')
        
        fig, ax = plt.subplots(1,1,figsize=(10,5))
        ds_his_sel.plot.line(ax=ax, x='time')        
        ax.legend(ds_his_sel[var_loc].to_series(),loc=1,fontsize=8) #optional, to change legend location        
        if i==0: 
            plt.ylim([-2,2])
        plt.savefig(r'figs_output/fig_his_'+var_dat[i]+'.png')
        plt.close()
        
    ''' _map.nc file  -> prepare'''
    # fn        = r'DATA/DFM_map.nc'
    # ds_map    = xr.open_dataset(fn)
    # var       = [v for v in list(ds_map.variables.keys()) if 'Mesh2d' in v]
    # ds_map    = ds_map[var]
    # ds_map.to_netcdf('DATA/DFM_map_2d.nc')
    # print(var)
    
    ''' _map.nc file  -> grid'''
    crs       = 'EPSG:28992'
    line      = np.array([[ 140000, 443000 ],[145000, 443000]])
    fn        = r'DATA/DFM_map_2d.nc'
    ds_map    = xr.open_dataset(fn)
    var       = [v for v in list(ds_map.variables.keys()) if 'Mesh2d' in v]
    print(var)
    ds_map    = ds_map[var]
    ds_map    = xu.core.wrap.UgridDataset(ds_map)
    fig, ax   = plt.subplots(figsize=(10,4))
    pc        = ds_map.grids[0].plot(edgecolor='crimson', linewidth=0.5,alpha=0.5)
    # Esri.WorldImagery
    # CartoDB.Voyager
    ctx.add_basemap(ax=ax, source=ctx.providers.CartoDB.Voyager, crs=crs, attribution=False)
    ax.set_aspect('equal')
    ax.plot(line[:,0],line[:,1],'b',linewidth=3)
    plt.savefig(r'figs_output/fig_map_grid.png')
    plt.close()
    
    ''' _map.nc file  -> bedlevel'''
    crs       = 'EPSG:28992'
    fn        = r'DATA/DFM_map_2d.nc'
    var       = 'Mesh2d_flowelem_bl'
    raster_res= 500
    ds_map    = xr.open_dataset(fn)
    ds_map    = xu.core.wrap.UgridDataset(ds_map)
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,7),sharex=True,sharey=True)
    pc = ds_map[var].ugrid.plot(ax=ax1, cmap='jet')
    pc = ds_map[var].ugrid.plot.contourf(ax=ax2, levels=11, cmap='jet')
    pc = ds_map[var].ugrid.plot.contour(ax=ax3, levels=11, cmap='jet', add_colorbar=True)
    bl_raster = dfmt.rasterize_ugrid(ds_map[[var]],resolution=raster_res) #rasterize ugrid
    pc = bl_raster[var].plot(ax=ax4, cmap='jet') #plot with non-ugrid method
    tit= ['plot','contourf','contour','plot rasterized']
    i=0
    for ax in (ax1,ax2,ax3,ax4):
        ctx.add_basemap(ax=ax, source=ctx.providers.Esri.WorldImagery, crs=crs, attribution=False)    
        ax.set_title(tit[i])
        i=i+1
    plt.savefig(r'figs_output/fig_map_bedlevel.png')
    plt.close()
    
    ''' _map.nc file -> different crs'''
    crs       = 'EPSG:28992'
    fn        = r'DATA/DFM_map_2d.nc'
    to_crs    = 'EPSG:4326'
    ds_map    = xr.open_dataset(fn)
    ds_map    = xu.core.wrap.UgridDataset(ds_map)
    ds_map.ugrid.set_crs(crs)
    uds_map_conv = ds_map.ugrid.to_crs(to_crs)
    fig, ax = plt.subplots()
    uds_map_conv["Mesh2d_waterdepth"].isel(time=0).ugrid.plot(ax=ax)
    ctx.add_basemap(ax=ax, source=None, crs=to_crs, attribution=False)
    plt.savefig(r'figs_output/fig_map_diff_crs.png')
    plt.close()
    
    ''' _map.nc file -> slice via x/y'''
    fn          = r'DATA/DFM_map_2d.nc'
    ds_map      = xr.open_dataset(fn)
    ds_map      = xu.core.wrap.UgridDataset(ds_map)
    sel_slice_x1, sel_slice_y1 = slice(140000,145000), slice(440000,445000)
    sel_slice_x2, sel_slice_y2 = slice(140000,145000), slice(443000,443500)
    ds_map_sel1 = ds_map.ugrid.sel(x=sel_slice_x1,y=sel_slice_y1)
    ds_map_sel2 = ds_map.ugrid.sel(x=sel_slice_x2,y=sel_slice_y2)
    
    fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(15,10))
    ds_map_sel1['Mesh2d_flowelem_bl'].ugrid.plot(ax=ax1)
    ds_map_sel2['Mesh2d_flowelem_bl'].ugrid.plot(ax=ax2)
    ds_map_sel1['Mesh2d_flowelem_bl'].plot(ax=ax3)
    ds_map_sel2['Mesh2d_flowelem_bl'].plot(ax=ax4)
    ax5.plot(ds_map_sel1['Mesh2d_flowelem_bl'].Mesh2d_face_x.values,ds_map_sel1['Mesh2d_flowelem_bl'].Mesh2d_face_y.values,'.')
    ax6.plot(ds_map_sel2['Mesh2d_flowelem_bl'].Mesh2d_face_x.values,ds_map_sel2['Mesh2d_flowelem_bl'].Mesh2d_face_y.values,'.')
    ax1.set_title('Box slice')
    ax2.set_title('Line slice')
    for ax in (ax5,ax6):
        ax.set_xlabel('x coord')
        ax.set_ylabel('y coord')    
    plt.savefig(r'figs_output/fig_map_slice_xy.png')
    plt.close()
    
    '''_map.nc file -> water level filter dry cells'''
    crs         = 'EPSG:28992'
    fn          = r'DATA/DFM_map_2d.nc'
    var         = 'Mesh2d_s1' # water level
    ds_map      = xr.open_dataset(fn)
    ds_map      = xu.core.wrap.UgridDataset(ds_map)
    bool_drycells = ds_map[var]==ds_map['Mesh2d_flowelem_bl']
    ds_map['mesh2d_s1_filt'] = ds_map[var].where(~bool_drycells)
    
    #plot water level on map
    fig, ax = plt.subplots(figsize=(10,4))
    pc = ds_map['mesh2d_s1_filt'].isel(time=3).ugrid.plot(cmap='jet')
    ctx.add_basemap(ax=ax, source=ctx.providers.Esri.WorldImagery, crs=crs, attribution=False)
    plt.savefig(r'figs_output/fig_map_WL_excl_dry_cells.png')
    plt.close()
    
    '''_map.nc file -> velocity wet cells (magnitude and quiver)'''
    crs         = 'EPSG:28992'
    fn          = r'DATA/DFM_map_2d.nc'
    ds_map      = xr.open_dataset(fn)
    ds_map      = xu.core.wrap.UgridDataset(ds_map)
    raster_res  = 500
    scale       = 7.5
    
    uds_quiv           = ds_map.isel(time=200, mesh2d_nLayers=-2, nmesh2d_layer=-2, missing_dims='ignore')
    varn_ucx, varn_ucy = 'Mesh2d_ucx', 'Mesh2d_ucy'
    magn_attrs         = {'long_name':'velocity magnitude', 'units':'m/s'}
    uds_quiv['magn']   = np.sqrt(uds_quiv[varn_ucx]**2+uds_quiv[varn_ucy]**2).assign_attrs(magn_attrs)
    raster_quiv        = dfmt.rasterize_ugrid(uds_quiv[[varn_ucx,varn_ucy]], resolution=raster_res)
    
    #plot
    fig,ax = plt.subplots(figsize=(10,4))
    pc = uds_quiv['magn'].ugrid.plot(cmap='jet')
    raster_quiv.plot.quiver(x='Mesh2d_face_x',y='Mesh2d_face_y',u=varn_ucx,v=varn_ucy,color='w',scale=scale,add_guide=False)
    ctx.add_basemap(ax=ax, source=ctx.providers.Esri.WorldImagery, crs=crs, attribution=False)
    plt.xlim([115000,140000])
    plt.ylim([440000,465000])
    plt.savefig(r'figs_output/fig_map_velocity_quiver.png')
    plt.close()
    
    '''_map.nc file -> slice/sideview '''
    crs         = 'EPSG:28992'
    fn          = r'DATA/DFM_map_2d.nc'
    ds_map      = xr.open_dataset(fn)
    ds_map      = xu.core.wrap.UgridDataset(ds_map)
    line        = np.array([[ 140000, 443000 ],[145000, 443000]])
    uds_crs     = dfmt.polyline_mapslice(ds_map.isel(time=3), line)
    
    fig, ax = plt.subplots()
    uds_crs['Mesh2d_flowelem_bl'].ugrid.plot(cmap='jet')
    plt.savefig(r'figs_output/fig_map_side_view.png')
    plt.close()