# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:51:27 2024

@author: PetraH
"""

import os
# import  aspose.cells 
# from aspose.cells import Workbook
from htmlmerger import HtmlMerger
from pathlib import Path

class paths():
    root = os.getcwd()
    folder = os.path.join(root,'Output')
    
if __name__ == "__main__":
    figs = ['fig_WL_vs_streven','fig_WL_T1','fig_dQ_laterals','fig_Q_T1','fig_Q_stat_dQ']
    fn   = [os.path.join(paths.folder,f + '.html') for f in figs]
    
    # merger = HtmlMerger(files=Path("Output/").glob("*html"))  # result will be in ./merged.html
    merger = HtmlMerger(files=fn)  # result will be in ./merged.html
    merger.merge()