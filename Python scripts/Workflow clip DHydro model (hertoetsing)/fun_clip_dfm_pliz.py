# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Verwijder alle regels die horen bij bij hoogtelijnen met coordinaten die vallen buiten het interesse gebied (rectangular o.b.v. shapefile)
"""

import numpy as np
import os
from tqdm import tqdm

def clip_dfm_pliz(path,fn_in, area):
    print('clip ' + fn_in)
    
    # load data
    df_in  = open(os.path.join(path,fn_in))
    lines  = df_in.readlines()
    [xmin,ymin,xmax,ymax] = area.total_bounds
    
    # create new output file
    fn_out   = fn_in.replace('.pliz','_sel.pliz')
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    # start index blocks of info
    idx_blocks = [i for i in range(0,len(lines)) if '*hdsr_ke' in lines[i]]
    
    with tqdm(total=len(idx_blocks)-1,desc=fn_in) as pbar:
        for i in range(0,len(idx_blocks)):
            row_start = idx_blocks[i]
            if row_start != idx_blocks[-1]:
                row_end   = idx_blocks[i + 1 ]                 
            else:
                row_end = len(lines)
            
            idx_data = [i for i in np.arange(row_start,row_end) if len(lines[i].split(' '))>=10]
            XX = np.array([float(lines[i].split(' ')[0]) for i in idx_data])
            YY = np.array([float(lines[i].split(' ')[2]) for i in idx_data])
            cond   = (np.min(XX)>=xmin) & (np.max(XX)<=xmax) & (np.min(YY)>=ymin) & (np.max(YY)<=ymax)
                
            if cond:                         
                for ii in range(row_start,row_end):
                    df_out.write(lines[ii])
                            
            pbar.update(1)
    
    df_in.close()
    df_out.close()

    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
    
    