"""
Script om 2D grid en 1D2D links uit D-HYDRO model te verwijderen. Andere acties handmatig gedaan

Auteur:             Boyan Domhof (HDSR)
Laatste update:     11-09-2024
"""

import xarray as xr

branch_prefixes_to_keep = ['rijn_wl', 'noor_wl']

# load 1d2d network netcdf file
ds = xr.open_dataset('./dflowfm/network.nc')

# remove 2d grid
v = [v for v in ds.variables if v.lower().startswith('mesh2d')]
ds = ds.drop_vars(v)
v = [v for v in ds.dims if v.lower().startswith('mesh2d')]
ds = ds.drop_dims(v)

# remove 1d2d links
v = [v for v in ds.variables if v.lower().startswith('links')]
ds = ds.drop_vars(v)
v = [v for v in ds.dims if v.lower().startswith('links')]
ds = ds.drop_dims(v)

# save to new network.nc file
ds.to_netcdf('./dflowfm/network1d.nc')

x = 2
