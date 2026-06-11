# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Verwijder alle regels die horen bij bij:
    - kunstwerken die weggefilterd zijn (dus niet inbegrepen in structures.ini)
    - laterals voor knopen die niet inbegrepen zijn in bnd.ext of nodeFile.ini
     -randvoorwaarden voor knopen die niet inbegrepeen zijn in bnd.ext


Note: 
    In DFM_lateral_sources.bc staat de eenheid m3/s (met superscript). 
    Dat geeft een error als deze workflow in spyder gerund wordt (in clip_dfm_bc). In anaconda/miniforge gaat het goed.
    Dit wordt nu voorkomen in fun_clip_dfm_bc.py door de command ", errors='ignore'" tijdens het laden van bc bestanden.
    Hierdoor verandert de eenheid wel van [m3/s] naar [m/s].
    

"""

import numpy as np
import os
from tqdm import tqdm

    
def clip_dfm_bc(path,fn_in, groupname,file_ref):
    print('clip ' + fn_in)
    # load data
    try:
        df_in  = open(os.path.join(path,fn_in)) # soms crasht dit: bij DFM_lateral_sources.bc door superscript in eenheid m3/s
        lines  = df_in.readlines()
    except:
        df_in  = open(os.path.join(path,fn_in),encoding='utf-8', errors='ignore') # "fix" code crash: ignore error; only needed in spyder
        lines  = df_in.readlines()
        
    df_ref = open(os.path.join(path,file_ref[0]))
    lines_ref = df_ref.readlines()
    if fn_in == 'DFM_lateral_sources.bc':
        df_ref2 = open(os.path.join(path,file_ref[1]))
        lines_ref2 = df_ref2.readlines()
    
    # selection
    if fn_in == 'DFM_lateral_sources.bc':
        # file_ref1: 'nodeFile.ini'
        lines_ref = [l.replace('\t','') for l in lines_ref]
        l_sel     = [l for l in lines_ref if 'name=' in l]  
        id_sel    = np.array([l.split('=')[1].strip() for l in l_sel])
        
        # file_ref2: 'bnd.ext'
        lines_ref2 = [l.replace('\t','') for l in lines_ref2]
        l_sel2     = [l for l in lines_ref2 if 'id=' in l]   
        id_sel2    = np.array([l.split('=')[1].strip() for l in l_sel2])
    else:
        lines_ref = [l.replace('\t','') for l in lines_ref]
        if file_ref == 'bnd.ext':
            l_sel  = [l for l in lines_ref if 'nodeId=' in l]        
        elif file_ref == 'nodeFile.ini':
            l_sel  = [l for l in lines_ref if 'name=' in l]      
        else:
            l_sel  = [l for l in lines_ref if 'id                    = ' in l]                
        id_sel = np.array([l.split('=')[1].strip() for l in l_sel])
    
    
        
    # create new output file
    fn_out   = fn_in.replace('.bc','_sel.bc')
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    # start index blocks of info
    idx_blocks = np.where(np.array(lines)==groupname)[0] 
    
    # check for groupname types
    check = np.unique(np.array([l for l in lines if '[' in l]))
    # print(check)
    
    # keep first block
    for ii in range(0,idx_blocks[0]):
        df_out.write(lines[ii])
    
    with tqdm(total=len(idx_blocks)-1,desc=fn_in) as pbar:
        for i in range(0,len(idx_blocks)):
            row_start = idx_blocks[i]
            if row_start != idx_blocks[-1]:
                row_end   = idx_blocks[i + 1 ]
            else:
                row_end = len(lines)
            
            # keep next block under specific conditions
            categories = np.array([l.split('=')[0].replace(' ','').replace('\t','') for l in lines[row_start:row_end]])
            idx_crit   = row_start + np.where(categories=='name')[0][0] # row in block starting with "branchID"
            id_name    = lines[idx_crit].split('=')[1].replace(' ','')[:-1] # corresponding string
            cond       = id_name in id_sel
            
            if fn_in == 'DFM_lateral_sources.bc':
                cond       = (id_name in id_sel) | (id_name in id_sel2)
                
            if cond:                  
                for ii in range(row_start,row_end):
                    df_out.write(lines[ii])
                
            pbar.update(1)
    
    df_in.close()
    df_out.close()
    df_ref.close()
    
    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
    

    
