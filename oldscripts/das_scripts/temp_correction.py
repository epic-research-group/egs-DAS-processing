#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 07:22:30 2022

@author: spri902
"""

import os
import gc
import math
from nptdms import TdmsFile
import pandas as pd
import numpy as np

import time
import matplotlib.dates as mdates
from datetime import datetime as datetime
import matplotlib
%matplotlib
import matplotlib.pyplot as plt 
import pickle
from boreholes import parse_surf_boreholes

#%%
# glass coefficient of thermal expansion
cte = 5.5e-7
# change in strain rate for thermal expansion

# dEps/dt = alpha*(d/dt)deltaT

# dEps/dt      = strain rate
# alpha        = cte of slica glass
# d/dt(deltaT) = time rate of change of detlaT (DTS temps)
#%% read in DAS and DTS data
df_full = pd.read_pickle('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/maystim22_26_combined_full')
dts_list=[]
for path,subdir,files in sorted(os.walk('/home/spri902/EGS_Collab/4850/DTS/processed/data')):
    for file in sorted(files):
        with open(f'{file}','rb') as f:
            well_dict = pickle.load(f)
            dts_list.append(well_dict)
            
