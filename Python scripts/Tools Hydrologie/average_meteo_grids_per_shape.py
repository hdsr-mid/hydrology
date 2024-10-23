"""
Average meteo grid cell values per shape in given shapefile. Only ASCII grid files are supported

Author:         Boyan Domhof (HDSR)
Last update:    23/10/2024
"""

from glob import glob
from pathlib import Path
import rioxarray as rioxr
import xarray as xr
import rasterio
from rasterio.features import shapes
import geopandas as gpd
import pandas as pd
import numpy as np
from tqdm import tqdm
import tempfile
import shutil


# Settings
meteo_grid_dir = Path('H:\\DATA\\Service\\Boyan\\00_WIWB_Raster_data')
meteo_var = 'precipitation'
crs = 'EPSG:28992'
shapefile_path = Path('H:\\DATA\\Service\\Boyan\\02_Dynamisch_Peilbeheer_Polder_Zegveld\\1_GIS\\Afvoergebied_Zegveld.shp')
shapefile_id_field = 'NAAM'

# Do the thing!
# first get shapes from shapefile
print(f'Load shapefile {shapefile_path}')
print(f'Shape id field {shapefile_id_field}')
gdf = gpd.read_file(shapefile_path)

# initialize dataframe to store values in
times = {f: '' for f in glob(str(meteo_grid_dir / meteo_var) + '\\*') if f.lower().endswith('.asc')}
times = {f: pd.to_datetime(f.replace(meteo_var, '').split('_')[-1].split('.')[0]) for f in times.keys()}
df = pd.DataFrame(index=times.values(), columns=gdf[shapefile_id_field].values)

# calculate weights raster per shape for clipping raster values accounting for intersecting shape edges
print('Calculate raster weights per shape')
# make temp copy of original ascii file to avoid modifying it
ftemp = Path(tempfile.gettempdir()) / 'temp.asc'
shutil.copy2(list(times.keys())[0], ftemp)
raster = rasterio.open(ftemp, 'r+')
# set CRS on raster
raster.crs = rasterio.crs.CRS.from_string(crs)
# convert raster cells to polygons
gdf_raster = gpd.GeoDataFrame.from_features(
    [{'properties': {'raster_val': v}, 'geometry': s}
    for s, v in shapes(np.arange(np.prod(raster.read().shape)).reshape(raster.read().shape), transform=raster.transform)]
)
gdf_raster.set_crs(rasterio.crs.CRS.from_string(crs), inplace=True)
# calculate weights
for i, row in gdf.iterrows():
    gdf_raster[row[shapefile_id_field]] =  gdf_raster.intersection(gdf.geometry.values[0]).area / gdf_raster.geometry.area
# convert weights back to weighted raster per shape (as a 2D numpy array)
gdf['weights'] = None
for i, row in gdf.iterrows():
    gdf.loc[i, 'weights'] = gdf_raster[row[shapefile_id_field]].values.reshape(raster.read().shape)

# loop over each ascii file, get average file per shape and store in dataframe
for f, t in tqdm(times.items(), 
                 desc=f'Calculate average meteo value {meteo_var} for each shape in shapefile'):
    ds = rioxr.open_rasterio(f, masked=True, crs=crs).isel(band=0).squeeze()
    # get average value per shape
    for i, row in gdf.iterrows():
        # get proper average value by multiplying with the precalculated weights which account for smaller clipped cells on the shape edges
        m = row['weights']
        m[row['weights'] <= 0] = np.nan
        total = np.nansum(ds*m)
        area = np.nansum(~np.isnan(ds)*m)
        df.loc[t, row[shapefile_id_field]] = total / area if area > 0 else np.nan

# export dataframe with average values to csv file
fout = meteo_grid_dir / f'{meteo_var}_{df.index[0]:%Y%m%d}_{df.index[-1]:%Y%m%d}.csv'
print(f'Export average values to csv file {fout}')
# a bit of processing to clean the data 
df.sort_index(inplace=True)
for c in df.columns:
    df.loc[df[c] > 200, c] = np.nan
df = df.resample(df.index[1] - df.index[0]).asfreq()
# export
df.to_csv(fout, sep=';', index_label='datetime', na_rep=-999)

print('Done')
