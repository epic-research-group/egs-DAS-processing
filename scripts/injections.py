#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 04:37:12 2022

@author: spri902
"""
import os
import numpy
import matplotlib
import matplotlib.pyplot as plt 
import time
import matplotlib.dates as mdates
from datetime import datetime as datetime
import pickle
import pandas as pd
#%%
os.chdir('/home/spri902/EGS_Collab/4850/stimflow/')
injFiles = sorted(os.listdir('/home/spri902/EGS_Collab/4850/stimflow/'))
injDat = pd.concat((pd.read_csv(f,header=1,usecols=[0,17],\
                           parse_dates = [0],infer_datetime_format=True) \
                    for f in injFiles[0:3] if f.endswith('.csv')),axis=0)
injDat.rename(columns={'hh:mm:ss':'date','psig.3':'psig'},inplace=True)
tmp = pd.read_csv(injFiles[3],header=1,usecols=[0,16],\
                          parse_dates=[0],infer_datetime_format=True)
tmp.rename(columns={'hh:mm:ss':'date','psig.2':'psig'},inplace=True)
injDat = pd.concat((injDat,tmp),axis=0)
injDat.reset_index(drop = True,inplace = True)
injDat.set_index('date',inplace=True)
#%%
ax = injDat.loc[:,("psig")].plot()
ax.set_ylabel('Injection Pressure (psig)')

