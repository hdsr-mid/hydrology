# -*- coding: utf-8 -*-
"""
Last modified: 08-06-2026
author: Petra Hulsman

Doel: Clip tif bestand o.b.v. polygon
Vervang naam in initialFields.ini

"""

import os
import rasterio
from rasterio.mask import mask
from shapely.geometry import mapping

def get_updated_meta(f, out_img, out_transform):
    """ Raster meta-informatie"""
    out_meta = f.meta.copy()
    out_meta['transform'] = out_transform
    out_meta['compress'] = 'deflate'
    out_meta['height'] = out_img.shape[1]
    out_meta['width'] = out_img.shape[2]
    #out_meta['crs'] = {'init': 'epsg:28992'}
    out_meta['crs'] = None
    return out_meta

def clip_dfm_tif(path, area):
    # define paths
    fn_orig  = 'wd_0_v5_5m_all-10mNAP.tif'
    fn_new   = 'wd_clip.tif'
    path_orig= os.path.join(path,fn_orig)
    path_new = os.path.join(path,fn_new)
    
    # Clip the raster    
    tif        = rasterio.open(path_orig)
    geometries = area.geometry.apply(mapping)

    tif_clip, out_transform = mask(tif, geometries, crop=True)
    out_meta = get_updated_meta(tif, tif_clip, out_transform)            
    with rasterio.open(path_new, 'w', **out_meta) as dest:
        dest.write(tif_clip)
    tif.close()
        
    # rename file in 'initialFields.ini'
    fn_in    = 'initialFields.ini'
    df_in    = open(os.path.join(path,fn_in))
    lines    = df_in.readlines()

    # create new output file
    extension= fn_in.split('.')[1]
    fn_out   = fn_in.replace('.'+extension,'_sel.'+extension)
    path_out = os.path.join(path,fn_out)
    if os.path.exists(path_out): os.remove(path_out)
    df_out = open(path_out, 'w')
    
    for i in range(0,len(lines)):
        line = lines[i]
        if fn_orig in line:
            line = line.replace(fn_orig,fn_new)
        df_out.write(line)

    df_in.close()
    df_out.close()
    
    # remove old file, rename new
    os.remove(os.path.join(path,fn_in))
    os.remove(os.path.join(path,fn_orig))
    os.rename(os.path.join(path,fn_out),os.path.join(path,fn_in))
