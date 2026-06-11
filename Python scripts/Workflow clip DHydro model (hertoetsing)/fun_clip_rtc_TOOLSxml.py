# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Verwijder regeles die horen bij kunstwerken/observatiepunten die weggeknipt zijn

"""

import numpy as np
import os
from tqdm import tqdm


gebiedsregelingen_koppelingen = {'Gebiedsregeling_KrommeRijn': ['hdsr_ge_G3003_benstrooms', 'hdsr_st_ST3001_bovstrooms','hdsr_st_ST3003_bovstrooms','hdsr_st_ST2028_bovstrooms','RTC_amelisweerd_downstream',
                                                                'hdsr_ge_G3003','hdsr_st_ST3001','hdsr_st_ST3003','hdsr_st_ST2028','Inlaat_WijkbDuurstede'],
                     'Gebiedsregeling_Schalkwijk': ['hdsr_ge_G4026_benstrooms','hdsr_st_ST4950_bovstrooms','KW1026_downstream',
                                                    'hdsr_ge_G4026','hdsr_st_ST4950'],
                     'Gebiedsregeling_Utrecht_GHIJ': ['RTC_SL0029_upstream','RTC_SL0009_upstream','RTC_SL0015_upstream','hdsr_ge_G8000_bovstrooms','RTC_SL0029_downstream''RTC_Nieuwegracht_downstream','RTC_ADCPHooggraven_downstream',
                                                      'hdsr_sl_SL0029','hdsr_sl_SL0009','Kokers_OoginAl','hdsr_ge_G8000'],
                     'Regeling_Zuidersluis': ['RTC_SL0015_upstream','RTC_SL0016_upstream',
                                              'hdsr_sl_SL0016'],
                     'Gebiedsregeling_OudeRijn': ['RTC_SL0003_downstream','RTC_SL0008_upstream','hdsr_st_ST1341_bovstrooms','RTC_SL0005_upstream','KW4312_upstream','KW4370_upstream','RTC_I6048_upstream','RTC_I6043_upstream','RTC_I6101_upstream','RTC_Woerden_downstream','RTC_Teckop_downstream',
                                                  'hdsr_sl_SL0008','hdsr_sl_SL0005','hdsr_sl_SL0003','hdsr_st_ST1341','hdsr_af_I6098'],
                     'Sifons_ARK':['RTC_ARK_downstream','hdsr_du_SY0070','hdsr_du_SY0071','hdsr_du_SY6009']}

gebiedsregelingen_variabelen = {'Gebiedsregeling_KrommeRijn': [],
                     'Gebiedsregeling_Schalkwijk': ['peilverschil_biester_kerkeland'],
                     'Gebiedsregeling_Utrecht_GHIJ': ['hver_Utrecht','u_h3','u_h2','u_h1','calc_spuien'],
                     'Regeling_Zuidersluis': [],
                     'Gebiedsregeling_OudeRijn': ['h_verwacht','or_h6','or_h5','or_h4','or_h3','or_h2','or_h1'],
                     'Sifons_ARK':[],
                     }

gebiedsregelingen_overige = {'Gebiedsregeling_KrommeRijn': ['check_2_3_v2','check_2_3_v3a','check_2_3_v4','check_2_3_v3b','check_2_3_v5',
                                                               'regime2_G3003','regime3_G3003','regime2_ST3001','regime3_ST3001','regime2_ST3003','regime3_ST3003','regime2_ST2028','regime3_ST2028','regime2_inlaat','regime3_inlaat'],
                     'Gebiedsregeling_Schalkwijk': ['peilverschil_check','peilverschil_check2',
                                                    'interval_G4026', 'G4026_hoogwater','interval_ST4950_hoogwater','interval_ST4950'],
                     'Gebiedsregeling_Utrecht_GHIJ': ['Maalstop_SL0029','regime2_SL0029','regime3_SL0029','regime4_SL0029','regime2_SL0009','regime3_SL0009','regime4_SL0009','regime2_OoginAl','regime3_OoginAl','regime4_OoginAl','G8000_uit','regime2_G8000','regime3_G8000','regime4_G8000',
                                                      'check_regime2_OoginAl','check_regime3_OoginAl','check_regime2_SL0009','check_regime3_SL0009','check_maalstpo_SL0029','check_regime2_SL0029','check_regime3_SL0029','check_pomp','check_regime2_G8000','check_regime3_G8000'],
                     'Regeling_Zuidersluis': ['Spuien','Standaard_regeling (dicht)',
                                              'Check_start_spuien','check_stop_spuien','check_spuien_bovenpeil'],
                     'Gebiedsregeling_OudeRijn': ['Spuien_sluis_check1','Spuien_sluis_check2','Spuien_sluis_check3','Spuien_sluis_check4','Spuien_sluis_check5',
                                                  'check_regime2_1','check_regime2_2','check_regime2_3','check_regime3_3','check_regime2_4','check_regime3_4','check_regime2_5'],
                     'Sifons_ARK': ['SY0070_dicht','SY0072_dicht','SY6009_dicht','SY0070_open','SY0072_open','SY6009_open',
                                  'check_waterstand_ARK1','check_waterstand_ARK2','check_waterstand_ARK3']
                     }

def get_names(path, fn_in, groupname,identifyer):
    
    # load data
    path      = path.replace('rtc','dflowfm')
    df_in     = open(os.path.join(path,fn_in))
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


def conditioning(group_lines,id_keep,id_remove):
    check_keep1       = [x for x in id_keep['part1'] for l in group_lines if x in l]
    check_remove1    = [x for x in id_remove['part1'] for l in group_lines if x in l]
    
    if (len(check_keep1)>0) | (len(check_remove1)>0):    
        if (len(check_keep1)>0) & (len(check_remove1)==0):
            cond = 'True'
               
        elif (len(check_keep1)>0) & (len(check_remove1)>0):
            cond = 'Maybe'
        
        else:
            cond = 'False'
    else:       
        # check pas o.b.v. kunstwerk-code alleen (part2) als hierboven niets gevonden is. 
        # Dus kijk pas naar "G8000" als niets gevonden is voor "hdsr_ge_G8000" (-> voor gevallen als check_regime2_G8000)
        check_keep2      = [x for x in id_keep['part2'] for l in group_lines if x in l]
        check_remove2    = [x for x in id_remove['part2'] for l in group_lines if x in l]
        if (len(check_keep2)>0) & (len(check_remove2)==0):
            cond = 'True'
               
        elif (len(check_keep2)>0) & (len(check_remove2)>0):
            cond = 'Maybe'
        
        else:
            cond = 'False'
        
    return cond

def clip_rtc_TOOLSxml(path,path_input,fn_in,group_start,group_end):
    print('clip ' + fn_in)
    
    # load data
    df_in    = open(os.path.join(path,fn_in))
    lines    = df_in.readlines()
    
    # find id of structures/observation points that are excluded
    strucs1   = get_names(path=path, fn_in = 'structures.ini', groupname = '[Structure]\n',identifyer='id')
    strucs2   = get_names(path=os.path.join(path_input,'rtc'), fn_in = 'structures.ini', groupname = '[Structure]\n',identifyer='id')
    strucs2   = [s for s in strucs2 if s not in strucs1] # alle kunstwerken die verwijderd zijn    
    obs1      = get_names(path=path, fn_in = 'obsFile1D_obs.ini', groupname = '[ObservationPoint]\n',identifyer='name')
    obs2      = get_names(path=os.path.join(path_input,'rtc'), fn_in = 'obsFile1D_obs.ini', groupname = '[ObservationPoint]\n',identifyer='name')
    obs2      = [s for s in obs2 if s not in obs1] # alle observationpoints die verwijderd zijn    
    obs1.remove('RTC_zomerwinter_downstream') # dit wordt overal gebruikt, dus gebruik de kunstwerkcode ipv dit observatiepunt
    
    strucs1b  = [i.split('_')[-1] for i in strucs1 if 'hdsr_' in i] # voor combinatie gevallen zoals check_regime2_G8000    
    strucs2b  = [i.split('_')[-1] for i in strucs2 if 'hdsr_' in i] # voor combinatie gevallen zoals check_regime2_G8000        
    
    # KrommeRijn: obsSL0015 included, strucs SL0015 excluded...
    id_remove = {'part1':strucs2 + obs2, 'part2':strucs2b}
    id_keep   = {'part1':strucs1 + obs1, 'part2':strucs1b}
    
    # add internally defined variables    
    for x in gebiedsregelingen_koppelingen.keys():
        check = [y for y in gebiedsregelingen_koppelingen[x] if y in id_keep['part1']]
        if len(check)>0:
            if len(check) == len(gebiedsregelingen_koppelingen[x]):
                id_keep['part1'] = id_keep['part1'] + gebiedsregelingen_koppelingen[x] + gebiedsregelingen_variabelen[x] + gebiedsregelingen_overige[x]         
            else:
                # if 1 one of the items is missing in this chain, remove internally defined variable
                id_remove['part1'] = id_remove['part1'] + gebiedsregelingen_variabelen[x] + gebiedsregelingen_overige[x]
        else:
            # if 1 one of the items is missing in this chain, remove internally defined variable
            id_remove['part1'] = id_remove['part1'] + gebiedsregelingen_variabelen[x] + gebiedsregelingen_overige[x]
                
    # create new output file
    extension= fn_in.split('.')[1]
    fn_out   = fn_in.replace('.'+extension,'_sel.'+extension)
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out   = open(path_out, 'w')
    
    # create new output file
    extension= fn_in.split('.')[1]
    fn_skip  = fn_in.replace('.'+extension,'_skipped.'+extension)
    path_skip= os.path.join(path,fn_skip)
    if os.path.exists(path_skip): os.remove(path_skip)
    df_skip  = open(path_skip, 'w')
    
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
            
            # save lines under specific conditions
            cond = conditioning(lines[row_start:row_end+1],id_keep,id_remove)
            if cond=='True':
                for ii in range(row_start,row_end+1):
                    df_out.write(lines[ii])
            elif cond=='Maybe':
                for ii in range(row_start,row_end+1):
                    df_skip.write(lines[ii])
            
            # check for intermediate data blocks not linked to structures                
            if (next_row_start != row_end+1) & (next_row_start != len(lines)):  
                for ii in range(row_end+1,next_row_start):
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
                    df_out.write(lines[ii])
            
            pbar.update(dx)
    
    df_in.close()
    df_out.close()  
    df_skip.close()
    
    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))