# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 26/02/2026
virtual environment used: geo-env
"""

import shutil
import os



class general():
    
    yref      = 2014
    season    = 'winter'
    refmodel  = r'D:\workdir\Model_1D_edit1_' + season
    bc_data   = r'D:\workdir\_prep_data'
    outmodel  = r'D:\workdir\Model_runs\Model_1D_edit1_' + season  
    outdir    = r'D:\workdir\Model_runs' 
    
def copy_folder(y):    
    s = general.refmodel
    d = general.outmodel + str(y)
    
    if not os.path.exists(d):
        shutil.copytree(s, d)
    else:
        print("Already exists")

def copy_lats_file(y):
    s     = os.path.join(general.bc_data,'bestanden_' + general.season,'dflowfm','DFM_lateral_sources_'+str(y)+'.bc')
    d     = os.path.join(general.outmodel + str(y), 'dflowfm','DFM_lateral_sources.bc')
    
    if os.path.exists(d):
        os.remove(d)
    shutil.copyfile(s, d)

def edit_file_rr(y):
    # Load data
    path_in  = os.path.join(general.refmodel,'rr','DELFT_3B.INI')
    data_in  = open(path_in, encoding='utf-8', errors='ignore')
    lines    = data_in.readlines()
    yref     = general.yref
    path_out = os.path.join(general.outmodel + str(y),'rr','DELFT_3B.INI')
    
    # edit file    
    data_out = open(path_out, 'w', encoding='utf-8')
    for jj in range(0,len(lines)):
        if general.season=='winter':            
            if "StartTime='"+str(yref)+"/10/15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref),str(y))
            if "EndTime='"+str(yref+1)+"/04/15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref+1),str(y+1))            
            
        if general.season=='zomer':
            if "StartTime='"+str(yref)+"/04/15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref),str(y))
            if "EndTime='"+str(yref)+"/10/15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref),str(y))
        
        data_out.write(lines[jj])
                    

    # close files
    data_in.close()
    data_in.close()
    data_out.close()

     
def edit_file_rtc(y):
    # Load data
    path_in  = os.path.join(general.refmodel,'rtc','rtcRuntimeConfig.xml')
    data_in  = open(path_in, encoding='utf-8', errors='ignore')
    lines    = data_in.readlines()
    yref     = general.yref
    path_out = os.path.join(general.outmodel + str(y),'rtc','rtcRuntimeConfig.xml')
    
    # edit file
    data_out = open(path_out, 'w', encoding='utf-8')
    for jj in range(0,len(lines)):
        if general.season=='winter':            
            if str(yref)+"-10-15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref),str(y))
            if str(yref+1)+"-04-15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref+1),str(y+1))            
            
        if general.season=='zomer':
            if str(yref)+"-04-15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref),str(y))
            if str(yref)+"-10-15" in lines[jj]:
                lines[jj] = lines[jj].replace(str(yref),str(y))

        data_out.write(lines[jj])
                    

    # close files
    data_in.close()
    data_out.close()    


def edit_file_mdu(y):
    # Load data
    path_in  = os.path.join(general.refmodel,'dflowfm','DFM.mdu')
    data_in  = open(path_in, encoding='utf-8', errors='ignore')
    lines    = data_in.readlines()
    yref     = general.yref
    path_out = os.path.join(general.outmodel + str(y),'dflowfm','DFM.mdu')
    
    # edit file
    data_out = open(path_out, 'w', encoding='utf-8')
    for jj in range(0,len(lines)):
        if "RefDate           = "+str(yref) in lines[jj]:
            lines[jj] = lines[jj].replace(str(yref),str(y))          
        
        data_out.write(lines[jj])
                    

    # close files
    data_in.close()
    data_out.close()  
    
def edit_file_bat(y):
    
    # edit combi file
    path_in  = os.path.join(general.outdir,'run'+str(general.yref)+'.bat')
    data_in  = open(path_in, encoding='utf-8', errors='ignore')
    lines    = data_in.readlines()
    yref     = general.yref
    path_out = path_in.replace(str(yref),str(y))
    
    # edit file
    data_out = open(path_out, 'w', encoding='utf-8')
    for jj in range(0,len(lines)):
        if str(yref) in lines[jj]:
            lines[jj] = lines[jj].replace(str(yref),str(y))          
        
        data_out.write(lines[jj])
                    

    # close files
    data_in.close()
    data_out.close()  
    
    
if __name__ == "__main__":
    
    for y in range(2018,2019):
        print(y)
        # create new folders
        copy_folder(y)
        
        # replace file
        copy_lats_file(y)
        
        # replace year in file
        edit_file_rr(y)
        edit_file_rtc(y)
        edit_file_mdu(y)
    
    
        edit_file_bat(y)
        
