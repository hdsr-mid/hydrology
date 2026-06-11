# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Clip rr/.bc op basis van id_keep
    

"""
import os
import numpy as np
from tqdm import tqdm

def clip_rr_bc(path, fn_in, groupname, id_keep)  :
    print('clip ' + fn_in)
    # load data
    df_in  = open(os.path.join(path,fn_in))
    lines  = df_in.readlines()
    
    # create new output file
    fn_out   = fn_in.replace('.bc','_sel.bc')
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    # start index blocks of info
    idx_blocks = np.where(np.array(lines)==groupname)[0] 
    
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
            cond       = id_name in id_keep
            
            if cond:                  
                for ii in range(row_start,row_end):
                    df_out.write(lines[ii])
                
            pbar.update(1)
    
    df_in.close()
    df_out.close()
    
    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))

    
