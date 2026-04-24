# -*- coding: utf-8 -*-
"""

@author     : MaaS-user, Petra Hulsman
Last update : 19/02/2026
virtual environment used: geo-env

"""

import datetime
import pandas as pd
import numpy as np

t0 = datetime.date(2011,1,1)
t1 = datetime.date(2024,1,1)

date_range = pd.date_range(t0,t1)

years   = [str(d.year).zfill(2) for d in date_range]
months  = [str(d.month).zfill(2) for d in date_range]
days    = [str(d.day).zfill(2) for d in date_range]
values  = ['0' for i in np.zeros(len(date_range))]

df = pd.DataFrame({'year':years,'months':months,'days':days,'value':values})

df.to_csv('D:/workingdir/3_Output/ModelData/evap.txt',index=False, sep=' ')