# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Inhoud bestand:
    - X - Y - Z waarden
    
Doel: Verwijder alle regels met coordinaten die vallen buiten het interesse gebied (rectangular o.b.v. shapefile):
    

"""

import os
from tqdm import tqdm

def clip_dfm_xyz(path,fn_in, area):
    print('clip ' + fn_in)
    
    # load data
    df_in  = open(os.path.join(path,fn_in))
    lines  = df_in.readlines()
    [xmin,ymin,xmax,ymax] = area.total_bounds
    
    # create new output file
    fn_out   = fn_in.replace('.xyz','_sel.xyz')
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    with tqdm(total=len(lines),desc=fn_in) as pbar:
        for i in range(0,len(lines)):
            # keep next block under specific conditions
            line   = lines[i]
            values = [l for l in line.strip().split(' ') if l!='']
            x      = float(values[0])
            y      = float(values[1])
            cond   = (x>=xmin) & (x<=xmax) & (y>=ymin) & (y<=ymax)
                
            if cond:                  
                df_out.write(lines[i])
                            
            pbar.update(1)
    
    df_in.close()
    df_out.close()

    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
    
    