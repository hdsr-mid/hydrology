# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Clip rr/.3b op basis van id_keep
    

"""
import os
from tqdm import tqdm


def clip_rr_3b(path,fn_in, id_keep):
    print('clip ' + fn_in)
    
    # load data
    df_in  = open(os.path.join(path,fn_in))
    lines  = df_in.readlines()
    
    # create new output file
    extension= fn_in.split('.')[1]
    fn_out   = fn_in.replace('.'+extension,'_sel.'+extension)
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    for i in tqdm(range(0,len(lines))):
        # keep line under specific conditions
        line    = lines[i]       
        id_name = line.split(' ')[2].replace("'",'')
        cond    = id_name in id_keep
            
        if cond:        
            df_out.write(lines[i])
    
    df_in.close()
    df_out.close()
    
    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
