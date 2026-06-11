# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Selecteer de knopen die (via-via) gelinked zijn aan dflowfm-knopen -> id_keep
Selecteer lijnen waarin deze knopen voorkomen
    

"""

import numpy as np
import os
from tqdm import tqdm

def fun_get_lat_nodes(path):
    
    # load data
    path      = path.replace('rr','dflowfm')
    fn_in     = 'bnd.ext'
    groupname = '[Lateral]\n'
    df_in     = open(os.path.join(path,fn_in))
    lines     = df_in.readlines()

    # start index blocks of info
    idx_blocks = np.where(np.array(lines)==groupname)[0] 

    id_names = []
    for i in range(0,len(idx_blocks)):
        row_start = idx_blocks[i]
        if row_start != idx_blocks[-1]:
            row_end   = idx_blocks[i + 1 ]                 
        else:
            row_end = len(lines)
        
        # keep next block under specific conditions
        categories = np.array([l.split('=')[0].replace(' ','').replace('\t','').replace('\n','') for l in lines[row_start:row_end]])
        idx_crit   = row_start + np.where(categories=='id')[0][0] # row in block starting with "branchID"
        id_name    = lines[idx_crit].split('=')[1].replace(' ','')[:-1] # corresponding string
        if 'AW' not in id_name: # deze zitten niet in 3B_LINK.TP
            id_names   = np.append(id_names,id_name)           

    df_in.close()
    return id_names


def clip_rr_3bLinks(path,fn_in):
    print('clip ' + fn_in)
    
    # load data
    df_in  = open(os.path.join(path,fn_in))
    lines  = df_in.readlines()
    
    # create new output file
    fn_out   = fn_in.replace('.TP','_sel.TP')
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    # get id-names to keep
    id_names1 = fun_get_lat_nodes(path) # get dflowfm node names (link to rr)
    id_names2 = ['pav_'+i.replace('_lat','') for i in id_names1 if i.startswith('node_')] # get rr node names (starting with dflowfm-nodes)
    id_keep   = np.append(id_names1,id_names2) 
    id_names3 = [i.replace('lat_','') for i in id_keep if 'wwtp' in i] # get connecting wwtp
    id_keep   = np.append(id_keep,id_names3) 
    
    for i in tqdm(range(0,len(lines))):
        # keep next block under specific conditions
        line   = lines[i]
        values = [l for l in line.strip().split(' ') if l!='']
        id_i1  = values[13].replace("'",'')
        id_i2  = values[15].replace("'",'')
        cond   = (id_i1 in id_keep) or (id_i2 in id_keep)
        
        if cond:                  
            df_out.write(lines[i])            
                
            if id_i2 not in id_keep:
                id_keep = np.append(id_keep,id_i2)
                
    
    df_in.close()
    df_out.close()

    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
    
    return id_keep
    
    