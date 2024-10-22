"""
Download meteo grids from WIWB (precipitation and Makkink evapotranspiration) for a given period in monthly increments

Author:         Boyan Domhof (HDSR)
Last update:    22/10/2024
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.append(r"Z:\\")
from hdsrhipy import Meteorology


# Settings
period = ('2010-01-01', '2024-01-01')
dt = 'H'
wiwb_credentials = ('api-wiwb-hdsr-hydrologen', 'IftWBsmfwleYrJ4RUWcjobzV3QazxoBz')
download_path = Path(r'Z:\\meteo')

# Do the thing!
meteo = Meteorology()
times = pd.date_range(period[0], period[1], freq=dt)
# request data in monthly increments
for y in np.unique(times.year):
    for m in np.unique(times[times.year == y].month):
        print(f'Downloading data for month {y:04d}-{m:02d}')
        # correct for the transition to new year
        y1, m1 = y, m
        y2 = y + 1 if (m + 1) > 12 else y
        m2 = m + 1 if (m + 1) <= 12 else 1
        if (y < 2019):
            # Meteobase neerslag t/m 2018
            meteo.download_from_WIWB(
                credentials=wiwb_credentials, 
                datasource='Meteobase.Precipitation', 
                variable='precipitation', 
                start=f'{y1:04d}{m1:02d}01000000', 
                end=f'{y2:04d}{m2:02d}01000000', 
                timestep=f'1{dt}', 
                download_path=download_path
            )
        else:
            # IRC vanaf 2019
            meteo.download_from_WIWB(
                credentials=wiwb_credentials, 
                datasource='Knmi.International.Radar.Composite.Final.Reanalysis', 
                variable='precipitation', 
                start=f'{y1:04d}{m1:02d}01000000', 
                end=f'{y2:04d}{m2:02d}01000000', 
                timestep=f'1{dt}', 
                download_path=download_path
            )
        # # Makkink verdamping
        # meteo.download_from_WIWB(
        #     credentials=wiwb_credentials, 
        #     datasource='Meteobase.Evaporation.Makkink', 
        #     variable='evaporation', 
        #     start=f'{y1:04d}{m1:02d}01000000', 
        #     end=f'{y2:04d}{m2:02d}01000000', 
        #     timestep='1D',  # only daily values available
        #     download_path=download_path
        # )