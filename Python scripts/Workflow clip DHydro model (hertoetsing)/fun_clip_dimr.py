# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Verwijder regels die horen bij:
    - laterals niet inbegrepen in dflowfm/bnd.ext
    - rtc knopen niet inbegrepen in rtcToolsConfig.xml
    
"""

import numpy as np
import os
from tqdm import tqdm

def get_names(path, fn_in, groupname,identifyer):
    
    # load data
    df_in     = open(os.path.join(path,'dflowfm',fn_in))
    lines     = df_in.readlines()

    # start index blocks of info
    idx_blocks = np.where(np.array(lines)==groupname)[0] 

    id_names = []
    with tqdm(total=len(idx_blocks)-1,desc=fn_in) as pbar:
        for i in range(0,len(idx_blocks)):
            row_start = idx_blocks[i]
            if row_start != idx_blocks[-1]:
                row_end   = idx_blocks[i + 1 ]                 
            else:
                row_end = len(lines)
            
            # keep next block under specific conditions
            categories = np.array([l.split('=')[0].replace(' ','').replace('\t','').replace('\n','') for l in lines[row_start:row_end]])
            idx_crit   = row_start + np.where(categories==identifyer)[0][0] # row in block starting with "branchID"
            id_name    = lines[idx_crit].split('=')[1].replace(' ','')[:-1] # corresponding string
            id_names   = np.append(id_names,id_name)           
            pbar.update(1)

    df_in.close()
    return list(id_names)

def get_RTCTools_names(path):
    
    # load data
    fn_in     = 'rtcToolsConfig.xml'
    df_in     = open(os.path.join(path,'rtc',fn_in))
    lines     = df_in.readlines()
    
    # select lines
    selection = [l for l in lines if '[' in l and ']' in l]
    id_keep   = ['[' + l.split('[')[1].split('<')[0].split('"')[0] for l in selection]
    return id_keep

def fun_clip_dimr(path,fn_in,group_start,group_end):
    print('clip ' + fn_in)
    
    # load data
    df_in   = open(os.path.join(path,fn_in))
    lines   = df_in.readlines()
    
    # get id's
    lats    = get_names(path,'bnd.ext','[Lateral]\n','id')
    rtc     = get_RTCTools_names(path) # RTC-links
    # strucs  = get_names(path,'structures.ini','[Structure]\n','id') 
    # obs     = get_names(path,'obsFile1D_obs.ini','[ObservationPoint]\n','name')
    id_keep = lats + rtc #+  obs  + strucs #
    
    # create new output file
    extension= fn_in.split('.')[1]
    fn_out   = fn_in.replace('.'+extension,'_sel.'+extension)
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    # save first rows
    row_start = 0
    row_end   = [j for j in range(row_start,len(lines)) if group_start in lines[j]][0]
        
    for ii in range(row_start,row_end):
        df_out.write(lines[ii])
    row_start = row_end
    
    # filter next rows
    with tqdm(total=len(lines),desc=fn_in) as pbar:
        while row_end < len(lines):
            # check for 1-liners            
            try:
                next_row_start   = [j for j in range(row_start+1,len(lines)) if group_start in lines[j]][0]                
            except:
                next_row_start   = len(lines)
            
            # find start/end row nr
            if next_row_start == row_start+1:
                row_end = row_start # 1-liner
            else:                
                row_end   = [j for j in range(row_start+1,len(lines)) if group_end in lines[j]][0]                
            # print(row_start,row_end,len(lines))  
            
            # keep next block under specific conditions
            id_check   = [x for x in id_keep for l in lines[row_start:row_end+1] if x in l]
            if len(id_check)>0:
                for ii in range(row_start,row_end+1):
                    df_out.write(lines[ii])
            
            # check for intermediate data blocks not linked to structures                
            if (next_row_start != row_end+1) & (next_row_start != len(lines)):                
                for ii in range(row_end+1,next_row_start):
                    df_out.write(lines[ii])
            
            # start of next block of data
            dx = next_row_start - row_start
            row_start = next_row_start
            
            # go through last couple of lines
            next_row_end     = [j for j in range(row_start+1,len(lines)) if group_end in lines[j]]            
            if len(next_row_end)==0:  
                row_start = row_end+1
                row_end   = len(lines)
                for ii in range(row_start,row_end):
                    if group_start in lines[ii]:
                        id_check   = [x for x in id_keep if x in lines[ii]]
                        if len(id_check)>0:
                            df_out.write(lines[ii])
                    else:
                        df_out.write(lines[ii])
                    
            
            pbar.update(dx)
    
    df_in.close()
    df_out.close()
    
    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))