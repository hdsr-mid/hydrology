# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: 
   Clip DHydro model voor een afvoergebied
    
    
Stappen:
- Clip branches-shapefile o.b.v. afvoergebied-polygoon (shapefile)
- Clip dflowfm/network.nc o.b.v. geclipte branches
- Clip bijbehorende kuntwerken, profiel info, laterals, observation points, boundary conditions, initial conditions
- Clip bedlevel.xyz en hoogtelijnen (.pliz) o.b.v. box-coordinaten afgeleid van het afvoergebied-polygoon (shapefile)
- Clip tif o.b.v. polygon
- Clip RR-knopen
- Clip RTC-knopen/sturing
- Verwijder overbodige koppelingen in dimr

    
Note: 
    In DFM_lateral_sources.bc staat de eenheid m3/s (met superscript). 
    Dat geeft een error als deze workflow in spyder gerund wordt (in clip_dfm_bc). In anaconda/miniforge gaat het goed.
    Dit wordt nu voorkomen in fun_clip_dfm_bc.py door de command ", errors='ignore'" tijdens het laden van bc bestanden.
    Hierdoor verandert de eenheid wel van [m3/s] naar [m/s].
    
"""

import os
import shutil
import geopandas as gdp
import pandas as pd


from fun_clip_dfm_netcdf import clip_dfm_netcdf
from fun_clip_dfm_init import clip_dfm_ini
from fun_clip_dfm_bc import clip_dfm_bc
from fun_clip_dfm_xyz import clip_dfm_xyz
from fun_clip_dfm_pliz import clip_dfm_pliz
from fun_clip_dfm_tif import clip_dfm_tif

from fun_clip_rr_3bLinks import clip_rr_3bLinks
from fun_clip_rr_netcdf import clip_rr_netcdf
from fun_clip_rr_bc import clip_rr_bc
from fun_clip_rr_3b import clip_rr_3b

from fun_clip_rtc_xml import clip_rtc_xml
from fun_clip_rtc_TOOLSxml import clip_rtc_TOOLSxml
from fun_clip_dimr import fun_clip_dimr


class general():
    root      = r'H:\DATA\Service\Petra\Projecten\2026 DHydro_model_knippen'
    ref_model = 'Model1D2D_edit'
    # ref_model = 'Model_1D'
    
    # definitie interesse gebied o.b.v. afvoergebied code
    clip_code      = 'AFVGEB0012' # DeTol
    naam           = 'DeTol'
    
    buffer = 1 # [m] kies hier hoeveel buffer je om je gebied wilt hebben

if __name__ == "__main__":
        # define paths
        m_input  = os.path.join(general.root, 'Input', general.ref_model)
        
        # select region
        afvoergebieden = gdp.read_file(os.path.join(general.root, 'Input', 'Afvoergebieden_30032026.gpkg'))        
        if len(general.clip_code)==1:
            area           = afvoergebieden.loc[afvoergebieden.CODE==general.clip_code]
        else:
            idx_code       = [i for i in range(0,len(afvoergebieden)) if afvoergebieden.CODE.values[i] in general.clip_code]
            area           = afvoergebieden.iloc[idx_code]
        # area           = gdp.read_file(os.path.join(general.root, 'Input', 'DeTol.gpkg')) # je kunt ook een shapefile/geopackage gebruiken..
        
        # add buffer and change crs if needed
        area['geometry'] = area.buffer(general.buffer)
        if area.crs != 'EPSG:28992': 
            area = area.to_crs('EPSG:28992')
        
        # select branches
        shp         = gdp.read_file(os.path.join(general.root, 'Input', 'branches.gpkg'))
        shp.crs     = 'EPSG:28992'
        shp1        = gdp.sjoin(shp, area, predicate = 'within')[shp.columns]
        idx         = [i for i in range(0,len(shp)) if 'hdsr' not in shp.iloc[i].Name]
        shp2        = shp.iloc[idx]
        shp_clip    = pd.concat([shp1,shp2])
        
        
        # create new folder
        m_output = os.path.join(general.root, 'Output', general.ref_model+'_clipped_'+general.naam)
        s = os.path.join(m_input)  # source directory
        d = os.path.join(m_output)  # destination directory
        shutil.copytree(s, d) # copy folder
        
        # save clipped network
        shp_clip.to_file(os.path.join(m_output,'branches_clip.gpkg'))
        
        print('--------------------------- clip model: dflowfm -----------------------------')
        
        path  = os.path.join(m_output,'dflowfm')
        
        clip_dfm_netcdf(path,shp_clip, area)
        
        clip_dfm_ini(path=path, fn_in='structures.ini'            , groupname='[Structure]\n')
        clip_dfm_ini(path=path, fn_in='obsFile1D_obs.ini'         , groupname='[ObservationPoint]\n')
        clip_dfm_ini(path=path, fn_in='crsloc.ini'                , groupname='[CrossSection]\n')
        clip_dfm_ini(path=path, fn_in='crsdef.ini'                , groupname='[Definition]\n')
        clip_dfm_ini(path=path, fn_in='Initialwaterlevels.ini'    , groupname='[Branch]\n')      
        clip_dfm_ini(path=path, fn_in='nodeFile.ini'              , groupname='[StorageNode]\n')
        clip_dfm_ini(path=path, fn_in='bnd.ext'                   , groupname=['[Boundary]\n','[Lateral]\n'])
        
        clip_dfm_bc(path=path, fn_in='FlowFM_structures.bc'       , groupname='[forcing]\n', file_ref=['structures.ini'])
        clip_dfm_bc(path=path, fn_in='DFM_lateral_sources.bc'     , groupname='[forcing]\n', file_ref=['nodeFile.ini','bnd.ext'])
        clip_dfm_bc(path=path, fn_in='DFM_boundaryconditions1d.bc', groupname='[forcing]\n', file_ref=['bnd.ext'])
        
        if '2D' in general.ref_model:
            clip_dfm_xyz(path=path, fn_in='bedlevel.xyz', area=area)
            clip_dfm_pliz(path=path, fn_in='objects.pli_fxw.pliz', area=area)        
            clip_dfm_tif(path, area)
        
        
        
        print('--------------------------- clip model: rr -----------------------------')
        path  = os.path.join(m_output,'rr')
        
        id_keep = clip_rr_3bLinks(path = path,fn_in = '3B_LINK.TP')
        
        clip_rr_netcdf(path = path, id_keep = id_keep)
        clip_rr_3b(path = path,fn_in = '3B_NOD.TP', id_keep = id_keep)    
        clip_rr_3b(path = path,fn_in = 'Paved.3b', id_keep = id_keep)    
        clip_rr_3b(path = path,fn_in = 'WWTP.3b', id_keep = id_keep)    
        clip_rr_bc(path=path, fn_in='BoundaryConditions.bc'       , groupname='[Boundary]\n', id_keep = id_keep)
        
        
        print('--------------------------- clip model: rtc -----------------------------')
        path  = os.path.join(m_output,'rtc')
        clip_rtc_TOOLSxml(path = path,path_input = m_input, fn_in = 'rtcToolsConfig.xml',group_start = ['    <rule>','    <trigger>'],group_end = ['    </rule>', '    </trigger>'])
        
        clip_rtc_xml(path = path, fn_in = 'rtcDataConfig.xml',group_start = ['<timeSeries id='],group_end = ['</timeSeries>'])
        clip_rtc_xml(path = path, fn_in = 'state_import.xml',group_start = ['<treeVectorLeaf id'],group_end = ['</treeVectorLeaf>'])
        clip_rtc_xml(path = path, fn_in = 'timeseries_import.xml',group_start = ['<series>'],group_end = ['</series>'])
        
        
        print('--------------------------- clip dimr -----------------------------')
        path  = m_output
        fun_clip_dimr(path = path,fn_in = 'dimr_config.xml',group_start = '<item>',group_end = '</item>')        
