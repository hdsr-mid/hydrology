# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 11:25:47 2025

@author: MaaS-user

Steps (when keeping 2D network and filtering 1D network):
    1. select network
       a. remove 2D components   
       b. select network_nEdges based on network_branch_id
       c. select network_nNodes based on network_edge_nodes 
       d. selecteer geo_nodes based on selected network_nEdges       
    2. select mesh1d
       a. selecteer mesh1d_nEdges based on mesh1d_edge_branch
       b. selecteer mesh1d_nNodes based on mesh1d_node_branch
    3. select links based on links_contact_id
    4. apply selection
        a. 1D network
        b. 2D components
    5. renumber indices to ensure it ranges from 0 to ...
	   a. network_edge_nodes
       b. mesh1d_edge_nodes
       c. mesh1d_edge_branch
       d. mesh1d_node_branch
       e. links
       f. links_contact_id
    6. save netcdf

"""

import numpy as np
import xarray as xr
import os
import geopandas as gdp
# np.set_printoptions(suppress=True)

def create_lookup(idx_sel):
    n_from       = np.unique(idx_sel)
    n_to         = np.arange(0,len(n_from))
    lookup_table = dict()
    for i in range(0,len(n_from)): lookup_table[n_from[i]] = n_to[i]        
    return lookup_table

def renumber(lookup_table,xds,var):
    data = xds[var].values
    
    if len(data.shape)==1:
        nr = np.array([lookup_table.get(x, x) for x in data])
        
    elif len(data.shape)==2:
        nr_1 = [lookup_table.get(x, x) for x in data[:,0]] 
        nr_2 = [lookup_table.get(x, x) for x in data[:,1]] 
        nr   = np.array([nr_1,nr_2]).T
    else:
        print('error renumbering')
    
    return nr
 
    
    
def clip_netcdf(path):
    print('clip netcdf')
    # load    
    xds  = xr.open_dataset(os.path.join(path,'network.nc'))
    shp  = gdp.read_file('shp/branches_sel.gpkg') 
        
    # 1. select network
    # a. remove 2D components
    var_sel = list(xds.head())
    if '_1D2D' not in path:   # hence, it's 1D     
        var_sel = [v for v in var_sel if '2d' not in v] # remove 2d components
        var_sel = [v for v in var_sel if 'link' not in v] # remove 1d2d links
    
    # b. remove 1D2Dlinks 1D2D model (some scenarios)
    if 'nolinks' in path:
        var_sel = [v for v in var_sel if 'link' not in v] # remove 1d2d links
    
    # b. selecteer network_nEdges o.b.v. network_branch_id (e.g. rivier afkorting inbegrepen)
    # content variable "network_branch_id"      : string, format: locatie_type_getal
    # content dimension "network_nEdges"        : numbers, network edge number
    
    # Option 1: based on name
    # nmax         = len(xds['network_branch_id'].values)
    # river_abbr   = ['ark', 'mark', 'rijn', 'noor']
    # idx_n_branch = np.array([i for i in range(0,nmax) if any(str(xds['network_branch_id'].values[i]).split('_')[0].replace("b'","") in x for x in river_abbr)])
    # print(len(idx_n_branch))
    
    # Option 2: based on shapefile
    # shapefile indexing: 0 - water authority channels, 1 - river (incl. connectors), 2 - tunn, 3 - unde, 4 - intersections
    nmax         = len(xds['network_branch_id'].values)
    if path.endswith('riv\dflowfm'):
        id_sel       = shp.loc[shp.Selection==1].Name.values 
    elif path.endswith('rivplus\dflowfm'):
        id_sel       = shp.loc[shp.Selection!=0].Name.values 
    else:
        print('keep entire 1D network')
        id_sel       = shp.Name.values        
    id_all       = [str(xds['network_branch_id'].values[i]).replace("b'","").replace("'","").strip() for i in range(0,nmax)]
    idx_n_branch = np.array([i for i in range(0,nmax) if id_all[i] in id_sel])
    # test         = np.unique(np.array([str(xds['network_branch_id'].values[i]).split('_')[0].replace("b'","") for i in idx_n_branch]))
    # print(test)
                    

    # c. select network_nNodes based on network_edge_nodes
    # content variable "network_edge_nodes"      : numbers, [nr_start_node, nr_end_node] -> begin/end node_nr of each network branch
    # content dimension "network_nNodes"         : numbers, network node number
    idx_n_nodes = np.unique(xds['network_edge_nodes'].values[idx_n_branch])
    
    # d. select geo_nodes based on selected network_nEdges
    # content variable "network_geom_node_count"    : numbers, nr of nodes in each branch
    # content variable "network_geom_x/y"           : numbers, coordinates of each geom node
    idx_geom = np.array([])
    for i in idx_n_branch:
        i0 = np.sum(xds['network_geom_node_count'].values[0:i])
        nr = xds['network_geom_node_count'].values[i]
        idx_geom = np.append(idx_geom,np.arange(i0,i0+nr)).astype(int)
    
    # 2. select mesh1d
    # a. select mesh1d_nEdges based on mesh1d_edge_branch
    # content variable "mesh1d_edge_branch"     : numbers, network branch number for each mesh1d edge (multiple mesh1d edges for each network branch)
    # content dimension "mesh1d_nEdges"         : numbers, mesh1d edge number
    branch_nr   = xds.mesh1d_edge_branch.values
    idx_m_edges = np.array([x for x in range(0,len(branch_nr)) if branch_nr[x] in idx_n_branch])
    
    # b. select mesh1d_nNodes based on mesh1d_node_branch
    # content variable "mesh1d_node_branch"     : numbers, network branch number for each mesh1d node (multiple mesh1d nodes for each network branch)
    # content dimension "mesh1d_nNodes"         : numbers, mesh1d node number
    branch_nr   = xds.mesh1d_node_branch.values
    idx_m_nodes = np.array([x for x in range(0,len(branch_nr)) if branch_nr[x] in idx_n_branch])
    
    if 'links' in var_sel:    
        # 3. select links based on links_contact_id
        # content variable "links"                      : numbers, [m1dnode_nr, pixel_nr]
        # content variable "links_contact_id"           : string,   m1dnode_nr_pixel_nr
        # content variable "links_contact_long_name"    : string,   m1dnode_nr_pixel_nr
        lmax      = len(xds.links_contact_id.values)
        mesh_link = xds.links.values[:,0]
        idx_links = np.array([i for i in range(0,lmax) if mesh_link[i] in idx_m_nodes])
    
    
        # 4. apply selection
        # a. 1D network
        xds_sel    = xds.isel(network_nEdges = idx_n_branch, network_nNodes = idx_n_nodes, 
                                    network_nGeometryNodes = idx_geom, 
                                    mesh1d_nEdges = idx_m_edges, mesh1d_nNodes = idx_m_nodes,
                                    links_nContacts = idx_links)
    else:
        # 4. apply selection
        # a. 1D network
        xds_sel    = xds.isel(network_nEdges = idx_n_branch, network_nNodes = idx_n_nodes, 
                                    network_nGeometryNodes = idx_geom, 
                                    mesh1d_nEdges = idx_m_edges, mesh1d_nNodes = idx_m_nodes)
    # 4. apply selection
    # b. 2D components
    xds_sel    = xds_sel[var_sel]
    
    # 5. renumber indices to ensure it ranges from 0 to ...
    # a. network_edge_nodes
    # content variable "network_edge_nodes"      : numbers, [nr_start_node, nr_end_node] -> begin/end node_nr of each network branch
    lookup_table           = create_lookup(idx_n_nodes)
    nr_network_edge_nodes  = renumber(lookup_table,xds_sel,'network_edge_nodes')
    
    # b. mesh1d_edge_nodes
    # content variable "mesh1d_edge_nodes"       : numbers, [nr_start_node, nr_end_node] -> begin/end node_nr of each mesh1d edge
    lookup_table          = create_lookup(idx_m_nodes)
    nr_mesh1d_edge_nodes  = renumber(lookup_table,xds_sel,'mesh1d_edge_nodes')
    
    # c. mesh1d_edge_branch
    # content variable "mesh1d_edge_branch"      : numbers, network branch number for each mesh1d edge (multiple mesh1d edges for each network branch)
    lookup_table          = create_lookup(idx_n_branch)
    nr_mesh1d_edge_branch = renumber(lookup_table,xds_sel,'mesh1d_edge_branch')
    
    # d. mesh1d_node_branch
    # content variable "mesh1d_node_branch"      : numbers, network branch number for each mesh1d node (multiple mesh1d nodes for each network branch)
    lookup_table          = create_lookup(idx_n_branch)
    nr_mesh1d_node_branch = renumber(lookup_table,xds_sel,'mesh1d_node_branch')
    
    if 'links' in var_sel:
        # e. links
        # content variable "links"                      : numbers, [m1dnode_nr, pixel_nr]
        lookup_table          = create_lookup(idx_m_nodes)
        nr_links_p1           = np.array([lookup_table.get(x, x) for x in xds_sel['links'].values[:,0]])
        nr_links_p2           = xds_sel['links'].values[:,1]
        nr_links              = np.array([nr_links_p1,nr_links_p2]).T
        
        # f. links_contact_id
        # content variable "links_contact_id"           : string,   m1dnode_nr_pixel_nr
        # content variable "links_contact_long_name"    : string,   m1dnode_nr_pixel_nr
        nr_links_str           = [str(int(nr_links_p1[i])) + '_' + str(int(nr_links_p2[i])) for i in range(0,len(nr_links_p1))]    
        nr_links_str_short     = np.array([n.ljust(40).encode("utf-8") for n in nr_links_str])
        nr_links_str_long      = np.array([n.ljust(80).encode("utf-8") for n in nr_links_str])
    
    # 6. apply renumbering
    xds_sel['network_edge_nodes']      = xds_sel['network_edge_nodes'].copy(data=nr_network_edge_nodes)
    xds_sel['mesh1d_edge_nodes']       = xds_sel['mesh1d_edge_nodes'].copy(data=nr_mesh1d_edge_nodes)
    xds_sel['mesh1d_edge_branch']      = xds_sel['mesh1d_edge_branch'].copy(data=nr_mesh1d_edge_branch)
    xds_sel['mesh1d_node_branch']      = xds_sel['mesh1d_node_branch'].copy(data=nr_mesh1d_node_branch)    
    if 'links' in var_sel:
        xds_sel['links']                   = xds_sel['links'].copy(data=nr_links)
        xds_sel['links_contact_id']        = xds_sel['links_contact_id'].copy(data=nr_links_str_short)
        xds_sel['links_contact_long_name'] = xds_sel['links_contact_long_name'].copy(data=nr_links_str_long)
                
    # 8. save netcdf
    xds_sel.to_netcdf(os.path.join(path,'network_sel.nc'))
    
    # 9. remove old file, rename new
    xds.close()    
    os.remove(os.path.join(path,'network.nc'))
    os.rename(os.path.join(path,'network_sel.nc'),os.path.join(path,'network.nc'))
    
    # # 10. remove other 2D related files
    # if '_1D2D' not in path:
    #     # os.remove(os.path.join(path,'bedlevel.xyz'))
    #     os.remove(os.path.join(path,'objects.pli_fxw.pliz'))
    #     os.remove(os.path.join(path,'roughness_sample.xyz'))
    #     os.remove(os.path.join(path,'wd_0_v10.tif'))