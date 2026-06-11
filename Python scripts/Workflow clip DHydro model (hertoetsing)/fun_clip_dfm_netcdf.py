# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Clip network.nc

Steps (when keeping 2D network and filtering 1D network):
    1. select grid
    2. select network         
       a. select network_nEdges based on network_branch_id
       b. select network_nNodes based on network_edge_nodes 
       c. selecteer geo_nodes based on selected network_nEdges       
    3. select mesh1d
       a. selecteer mesh1d_nEdges based on mesh1d_edge_branch
       b. selecteer mesh1d_nNodes based on mesh1d_node_branch
    4. select links based on links_contact_id
    5. apply selection
        a. 1D network
        b. 2D components
    6. renumber indices to ensure it ranges from 0 to ...
	   a. network_edge_nodes
       b. mesh1d_edge_nodes
       c. mesh1d_edge_branch
       d. mesh1d_node_branch
       e. links
       f. links_contact_id
    7. apply renumbering    
    8. save netcdf
    9. remove old file, rename new

"""

import numpy as np
import xarray as xr
import os

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
        
    elif (len(data.shape)==2) & (data.shape[1]==2):
        nr_1 = [lookup_table.get(x, x) for x in data[:,0]] 
        nr_2 = [lookup_table.get(x, x) for x in data[:,1]] 
        nr   = np.array([nr_1,nr_2]).T
    
    elif (len(data.shape)==2) & (data.shape[1]==4):
        nr_1 = [lookup_table.get(x, x) for x in data[:,0]] 
        nr_2 = [lookup_table.get(x, x) for x in data[:,1]] 
        nr_3 = [lookup_table.get(x, x) for x in data[:,2]] 
        nr_4 = [lookup_table.get(x, x) for x in data[:,3]] 
        nr   = np.array([nr_1,nr_2,nr_3,nr_4]).T
    
    else:
        print('error renumbering')
    
    return nr
 
    
def clip_dfm_netcdf(path,shp, area):
    print('clip netcdf')
    # load    
    xds  = xr.open_dataset(os.path.join(path,'network.nc'))
    
    # 1. select grid
    # a. select grid nodes based on coordinates
    if 'Mesh2d' in list(xds.head()):
        [xmin,ymin,xmax,ymax] = area.total_bounds
        idx_node_grid = np.where((xds.Mesh2d_node_x.values>=xmin) & (xds.Mesh2d_node_x.values<=xmax) &
                              (xds.Mesh2d_node_y.values>=ymin) & (xds.Mesh2d_node_y.values<=ymax))[0]
        
        # b. select grid faces & edges based on grid nodes
        idx_face_grid = np.where(np.isin(xds.Mesh2d_face_nodes.values[:,0],idx_node_grid) & 
                                 np.isin(xds.Mesh2d_face_nodes.values[:,1],idx_node_grid) & 
                                 np.isin(xds.Mesh2d_face_nodes.values[:,2],idx_node_grid) & 
                                 np.isin(xds.Mesh2d_face_nodes.values[:,3],idx_node_grid))[0]
        idx_edge_grid = np.where(np.isin(xds.Mesh2d_edge_nodes.values[:,0],idx_node_grid) & 
                                 np.isin(xds.Mesh2d_edge_nodes.values[:,1],idx_node_grid))[0]
    
    # 2. select network
    # a. select network_nEdges based on network_branch_id
    # content variable "network_branch_id"      : string, format: locatie_type_getal
    # content dimension "network_nEdges"        : numbers, network edge number
    nmax         = len(xds['network_branch_id'].values)
    id_sel       = shp.Name.values  
    id_all       = [str(xds['network_branch_id'].values[i]).replace("b'","").replace("'","").strip() for i in range(0,nmax)]
    idx_n_branch = np.array([i for i in range(0,nmax) if id_all[i] in id_sel])
     
    # b. select network_nNodes based on network_edge_nodes
    # content variable "network_edge_nodes"      : numbers, [nr_start_node, nr_end_node] -> begin/end node_nr of each network branch
    # content dimension "network_nNodes"         : numbers, network node number
    idx_n_nodes = np.unique(xds['network_edge_nodes'].values[idx_n_branch])
    
    # c. select geo_nodes based on selected network_nEdges
    # content variable "network_geom_node_count"    : numbers, nr of nodes in each branch
    # content variable "network_geom_x/y"           : numbers, coordinates of each geom node
    idx_geom = np.array([])
    for i in idx_n_branch:
        i0 = np.sum(xds['network_geom_node_count'].values[0:i])
        nr = xds['network_geom_node_count'].values[i]
        idx_geom = np.append(idx_geom,np.arange(i0,i0+nr)).astype(int)
    
    # 3. select mesh1d
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
    
    if 'Mesh2d' in list(xds.head()):
        # 4. select links based on links_contact_id
        # content variable "links"                      : numbers, [m1dnode_nr, pixel_nr]
        # content variable "links_contact_id"           : string,   m1dnode_nr_pixel_nr
        # content variable "links_contact_long_name"    : string,   m1dnode_nr_pixel_nr
        lmax = len(xds.links_contact_id.values)
        ID1  = [int(x.decode("utf-8").strip().split('_')[0]) for x in xds['links_contact_id'].values]
        ID2  = [int(x.decode("utf-8").strip().split('_')[1]) for x in xds['links_contact_id'].values]
        idx_links = [i for i in range(0,len(ID1)) if ID1[i] in idx_m_nodes and ID2[i] in idx_face_grid]
    
    # 5. apply selection
    if 'Mesh2d' not in list(xds.head()):
        xds_sel    = xds.isel(network_nEdges = idx_n_branch, 
                              network_nNodes = idx_n_nodes, 
                              network_nGeometryNodes = idx_geom, 
                              mesh1d_nEdges = idx_m_edges, 
                              mesh1d_nNodes = idx_m_nodes,                               
                                    )
    else:
        xds_sel    = xds.isel(network_nEdges = idx_n_branch, 
                              network_nNodes = idx_n_nodes, 
                              network_nGeometryNodes = idx_geom, 
                              mesh1d_nEdges = idx_m_edges, 
                              mesh1d_nNodes = idx_m_nodes, 
                              links_nContacts = idx_links,
                              Mesh2d_nNodes = idx_node_grid,
                              Mesh2d_nFaces = idx_face_grid,
                              Mesh2d_nEdges = idx_edge_grid
                              )
        
    # 6. renumber indices to ensure it ranges from 0 to ...
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
    
    if 'Mesh2d' in list(xds.head()):
        # e. Mesh2d_edge_nodes
        # content variable "Mesh2d_edge_nodes"      : [nr_start_node, nr_end_node]
        lookup_table           = create_lookup(idx_node_grid)
        nr_2D_edge_nodes       = renumber(lookup_table,xds_sel,'Mesh2d_edge_nodes')
    
        # f. Mesh2d_face_nodes
        # content variable "Mesh2d_face_nodes"      : [nr_node at 4 corners of grid]
        lookup_table           = create_lookup(idx_node_grid)
        nr_2D_face_nodes       = renumber(lookup_table,xds_sel,'Mesh2d_face_nodes')
    
        # g. links
        # content variable "links"                      : numbers, [m1dnode_nr, pixel_nr]
        lookup_table1         = create_lookup(idx_m_nodes)
        lookup_table2         = create_lookup(idx_face_grid)
        nr_links_p1           = np.array([lookup_table1.get(x, x) for x in xds_sel['links'].values[:,0]])
        nr_links_p2           = np.array([lookup_table2.get(x, x) for x in xds_sel['links'].values[:,1]])
        nr_links              = np.array([nr_links_p1,nr_links_p2]).T
        
        # h. links_contact_id
        # content variable "links_contact_id"           : string,   m1dnode_nr_pixel_nr
        # content variable "links_contact_long_name"    : string,   m1dnode_nr_pixel_nr
        nr_links_str           = [str(int(nr_links_p1[i])) + '_' + str(int(nr_links_p2[i])) for i in range(0,len(nr_links_p1))]    
        nr_links_str_short     = np.array([n.ljust(40).encode("utf-8") for n in nr_links_str])
        nr_links_str_long      = np.array([n.ljust(80).encode("utf-8") for n in nr_links_str])
    
    # 7. apply renumbering
    xds_sel['network_edge_nodes']      = xds_sel['network_edge_nodes'].copy(data=nr_network_edge_nodes)
    xds_sel['mesh1d_edge_nodes']       = xds_sel['mesh1d_edge_nodes'].copy(data=nr_mesh1d_edge_nodes)
    xds_sel['mesh1d_edge_branch']      = xds_sel['mesh1d_edge_branch'].copy(data=nr_mesh1d_edge_branch)
    xds_sel['mesh1d_node_branch']      = xds_sel['mesh1d_node_branch'].copy(data=nr_mesh1d_node_branch)    
    if 'Mesh2d' in list(xds.head()):
        xds_sel['links']                   = xds_sel['links'].copy(data=nr_links)
        xds_sel['links_contact_id']        = xds_sel['links_contact_id'].copy(data=nr_links_str_short)
        xds_sel['links_contact_long_name'] = xds_sel['links_contact_long_name'].copy(data=nr_links_str_long)
        xds_sel['Mesh2d_edge_nodes']       = xds_sel['Mesh2d_edge_nodes'].copy(data=nr_2D_edge_nodes)
        xds_sel['Mesh2d_face_nodes']       = xds_sel['Mesh2d_face_nodes'].copy(data=nr_2D_face_nodes)
    
    # 8. save netcdf
    xds_sel.to_netcdf(os.path.join(path,'network_sel.nc'))
    
    # 9. remove old file, rename new
    xds.close()    
    os.remove(os.path.join(path,'network.nc'))
    os.rename(os.path.join(path,'network_sel.nc'),os.path.join(path,'network.nc'))
