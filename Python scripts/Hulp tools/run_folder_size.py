# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 08:43:57 2023

@author: PetraH
"""
import os
from tqdm import tqdm

if __name__ == "__main__":
    root='H:/DATA/Service/...' # paste path here
    size=0
    
    # get size
    for path, dirs, files in os.walk(root):
        for f in tqdm(files):
            fp = os.path.join(path, f)
            size += os.path.getsize(fp) * 10**(-9)
     
    # display size
    print('')
    print("Folder size: " + str(size) + ' GB')