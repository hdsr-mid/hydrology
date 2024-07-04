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
    figs = ['fig_WL_streefpeil','fig_WL_T1','fig_laterals_gem','fig_Q_T1','fig_Q_gem_dQ','fig_laterals_sat_gem','fig_Q_dir']
    fn   = [os.path.join(paths.folder,f + '.html') for f in figs]
    
    # merger = HtmlMerger(files=Path("Output/").glob("*html"))  # result will be in ./merged.html
    merger = HtmlMerger(files=fn)  # result will be in ./merged.html
    merger.merge()
    
    
    hti = Html2Image(browser='edge',size=(1500, 4500))
    with open('merged.html') as f:
        hti.screenshot(f.read(), save_as='merged.png')