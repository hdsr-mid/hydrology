# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:48:38 2024

@author: PetraH
"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

import io
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from urllib.request import urlopen, Request
from PIL import Image


root   = os.getcwd()


def image_spoof(self, tile): # this function pretends not to be a Python script
    url = self._image_url(tile) # get the url of the street map API
    req = Request(url) # start request
    req.add_header('User-agent','Anaconda 3') # add user agent to request
    fh = urlopen(req) 
    im_data = io.BytesIO(fh.read()) # get image
    fh.close() # close url
    img = Image.open(im_data) # open image with PIL
    img = img.convert(self.desired_tile_form) # set image format
    return img, self.tileextent(tile), 'lower' # reformat for cartopy



if __name__ == "__main__":

    band         = 'Mesh2d_waterdepth'
    data_var     = 'Waterdiepte (max) [m]'
    vmins        = 0
    vmaxs        = 2.5
    
    # Get data
    t          = 0
    xds        = xr.open_dataset(os.path.join(root,'data.nc'))
    time       = xds.time.values
    lon        = xds['Mesh2d_face_x_bnd'].values[:,0]
    lat        = xds['Mesh2d_face_y_bnd'].values[:,0]                          
    lon_i      = np.array([np.where(np.round(l,3)==np.unique(np.round(lon,3)))[0][0] for l in lon])
    lat_i      = np.array([np.where(np.round(l,3)==np.unique(np.round(lat,3))[::-1])[0][0] for l in lat])
    rmax       = np.max(lat_i)+1
    cmax       = np.max(lon_i)+1
    data       = xds[band].values
    
    # Convert data        
    data_t  = data[t]
    data_2d = np.zeros((rmax,cmax))*np.nan 
    data_2d[lat_i,lon_i] = data_t           
    data_2d[data_2d==0] = np.nan
    xds = xr.Dataset(data_vars=dict(Waterdiepte=(["y", "x"], data_2d)),coords=dict(x= np.unique(np.round(lon,3)),y=np.unique(np.round(lat,3))[::-1]))
    xds.rio.write_crs("epsg:28992", inplace=True)
    
    # Plot
    plt.figure(figsize=(12,7))
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.95)                    
    ax = plt.axes(projection=ccrs.epsg(28992))
    
    request = cimgt.OSM()
    ax.add_image(request, 10)    # 5 = zoom level
    xds['Waterdiepte'].plot.imshow('x','y', ax=ax, vmin = vmins, vmax = vmaxs, cmap = 'Blues', zorder=10, alpha=0.95)        
    plt.title(pd.to_datetime(time[t]) - pd.to_datetime(time[0]))
    
    plt.savefig('fig.png', dpi=300)
    plt.close()
    
    
    