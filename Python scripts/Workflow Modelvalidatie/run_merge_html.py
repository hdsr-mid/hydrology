# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:51:27 2024

@author: PetraH

environment: env_merge_html
"""

import os
import numpy as np
from htmlmerger import HtmlMerger
from pathlib import Path
from html2image import Html2Image

class paths():
    root = os.getcwd()
    folder = os.path.join(root,'Output')
    
if __name__ == "__main__":
    
    for i in range(2,3):
        figs = [f for f in os.listdir('Output') if f.startswith(str(i) + '_')]
        fn   = [os.path.join(paths.folder,f) for f in figs]
        
        # merger = HtmlMerger(files=Path("Output/").glob("*html"))  # result will be in ./merged.html
        merger = HtmlMerger(files=fn)  # result will be in ./merged.html
        merger.merge()
        
        
        hti = Html2Image(browser='edge',size=(1500, 4500))
        with open('merged.html') as f:
            hti.screenshot(f.read(), save_as=str(i) + '_merged_figures.png')
        
        os.rename('merged.html',str(i) + '_merged_figures.html')