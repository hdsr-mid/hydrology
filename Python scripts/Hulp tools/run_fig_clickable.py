# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 16:19:41 2024

@author: petra
"""

import numpy as np
from bokeh.plotting import figure, save
from bokeh.models import CustomJS, ColumnDataSource, TapTool
from bokeh.layouts import layout

if __name__ == "__main__":
    # Generate random data
    xx = np.linspace(0,10,500)
    y = np.sin(xx)
    barelements = [2,5,7]
    
    #create a dictionary for each bar element based the dummy data/equations you provided    
    #idea is that this dictionary will get passed to CustomJS and then used to replace the data in the datsource that's plotting the lines    
    src_dict = {'Line 1':{'x':xx,'y':y}
                ,'Line 2':{'x':xx*3,'y':y**2}
                ,'Line 3':{'x':xx/2,'y':y**3}}
    
    #make a columndatasource for your bar
    bar_src = ColumnDataSource({'x':list(src_dict.keys()),'y':barelements,'c':['red','blue','green']})
    
    #start with initial state for columndatasource, Line 1's keys
    src=ColumnDataSource(data=src_dict['Line 1'])
    
    figure_barchart = figure(height=250, title="Click a Bar",x_range=list(src_dict.keys()))
    bar_rend = figure_barchart.vbar(x='x', top='y', width=0.5, bottom=0, color='c',source=bar_src)
    
    p1 = figure(title='Line 1', height=150)
    
    line_rend = p1.line(x='x', y='y',line_color='red',source=src)
    
    #create the TapTool, add it to the bar figure, and enable it on initialization
    tap = TapTool(renderers=[bar_rend])
    figure_barchart.add_tools(tap)
    figure_barchart.toolbar.active_tap = tap
    
    cb = CustomJS(args=dict(bar_src=bar_src,src_dict=src_dict,src=src,p1=p1,line_rend=line_rend)
                  ,code='''
                  //get the index of the selected bar
                  //getting 0th index cuz we just want on index value
                  var sel_bar_i = bar_src.selected.indices[0]
                  // what is the line name associated with that index?
                  var line_id = bar_src.data['x'][sel_bar_i]
                  //now that we know that line id, we can update the line's datasource via the src_dict
                  src.data = src_dict[line_id]
                  //update the line renderer's color to match just for demo purposes
                  line_rend.glyph.line_color = bar_src.data['c'][sel_bar_i]
                  //update the line plot title too just for demo purposes
                  p1.title.text=line_id
                  src.change.emit()
                  p1.change.emit()
                  ''')
    #apply this callback to trigger whenever the indices of the selected glyph change
    bar_src.selected.js_on_change('indices',cb)             
    
    save(layout([figure_barchart,p1]),'fig_clickable.html')