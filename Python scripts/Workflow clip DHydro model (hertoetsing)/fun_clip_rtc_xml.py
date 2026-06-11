# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Verwijder in regel-blokken voor sturingen die niet inbegrepen zijn in rtcToolsConfig.xml

"""

import os
from tqdm import tqdm

def get_RTCTools_names(path):
    
    # load data
    fn_in     = 'rtcToolsConfig.xml'
    df_in     = open(os.path.join(path,fn_in))
    lines     = df_in.readlines()
    
    # select lines
    selection = [l for l in lines if '[' in l and ']' in l]
    id_keep   = ['[' + l.split('[')[1].split('<')[0].split('"')[0] for l in selection]
    return id_keep


def row_identifyer(group_start,line):
    if len(group_start[0])==len(group_start[0].strip()):
        if len(group_start)==1:
            cond = group_start[0] in line
        else:
            cond0 = group_start[0] in line
            cond1 = group_start[1] in line
            cond  = cond0 or cond1
    else:
        if len(group_start)==1:
            cond = line.startswith(group_start[0])
        else:
            cond0 = line.startswith(group_start[0])
            cond1 = line.startswith(group_start[1])
            cond  = cond0 or cond1
    return cond

def clip_rtc_xml(path,fn_in,group_start,group_end):
    print('clip ' + fn_in)
    
    # load data
    df_in    = open(os.path.join(path,fn_in))
    lines    = df_in.readlines()
    id_keep  = get_RTCTools_names(path)
        
    # create new output file
    extension= fn_in.split('.')[1]
    fn_out   = fn_in.replace('.'+extension,'_sel.'+extension)
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    # save first rows
    row_start = 0
    row_end   = [j for j in range(row_start,len(lines)) if row_identifyer(group_start,lines[j])][0]
        
    for ii in range(row_start,row_end):
        df_out.write(lines[ii])
    row_start = row_end
    
    # filter next rows
    with tqdm(total=len(lines),desc=fn_in) as pbar:
        while row_end < len(lines):
            # check for 1-liners            
            try:
                next_row_start   = [j for j in range(row_start+1,len(lines)) if row_identifyer(group_start,lines[j])][0]                
            except:
                next_row_start   = len(lines)
            
            # find start/end row nr
            if next_row_start == row_start+1:
                row_end = row_start # 1-liner
            else:                
                row_end   = [j for j in range(row_start+1,len(lines)) if row_identifyer(group_end,lines[j])][0]                            
            
            # keep next block under specific conditions
            id_check     = [x for x in id_keep for l in lines[row_start:row_end+1] if x in l]
            
            if len(id_check)>0:
                for ii in range(row_start,row_end+1):
                    df_out.write(lines[ii])
                    
            # check for intermediate data blocks not linked to structures                
            if (next_row_start != row_end+1) & (next_row_start != len(lines)):                 
                for ii in range(row_end+1,next_row_start):
                    # print('part2: ',lines[ii])
                    df_out.write(lines[ii])
            
            # start of next block of data
            dx = next_row_start - row_start
            row_start = next_row_start
            
            # go through last couple of lines
            next_row_end     = [j for j in range(row_start+1,len(lines)) if row_identifyer(group_end,lines[j])]            
            if len(next_row_end)==0:  
                row_start = row_end+1
                row_end   = len(lines)
                for ii in range(row_start,row_end):
                    if row_identifyer(group_start,lines[ii]):
                        id_check     = [x for x in id_keep if x in lines[ii] ]
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