#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 12:49:29 2022
make source receiver gathers for may 2018 stimulation
at EGS Collab experiment 1

@author: spri902
"""
import matplotlib
%matplotlib
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
import os
from obspy.io.seg2.seg2 import _read_seg2
import numpy as np
import pandas as pd
from phasepickers.parker_cfs import CASSM_chanDecode
from datetime import datetime as datetime
import multiprocessing as mp

#%%

directories = sorted(os.walk('/shared/mayCASSM/'))
dirpath='/shared/mayCASSM/'
dates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directories[0][1])]
siglen = 2000
OThydata = np.empty((len(dates),siglen,11))
cnt=0
postfix='dat'

#%%
def build_singleSRCgather(filename):
    global cnt
    global OThydata
    
    while cnt < len(dates):
        src=[0]
        recChan=slice(12,23)

        low=3000
        high=8000
        sthead = _read_seg2(filename, headonly=True)
        times = sthead[0].times()
        srcCh = CASSM_chanDecode(t=times, enc=sthead[93].data)
        if srcCh in src:
            dt = sthead[0].stats.delta
            fs = np.floor(sthead[0].stats.sampling_rate)
            st = sthead[0].stats.starttime
            et = sthead[0].stats.endtime
            sthead.trim(st, st + (dt*1999))
            OThy = sthead[recChan]
            OThy = OThy.filter('bandpass', freqmin=low,freqmax=high,
                           corners=4, zerophase=True)
        
            OThydata[cnt,:,:] = np.stack(tr.data for tr in OThy).T
        
            cnt += 1
        else:
            continue
    return OThydata


num_proc = 8
pool = mp.Pool(processes = num_proc)
proc = [pool.apply_async(build_singleSRCgather,args=[os.path.join(root,fn)]) 
        for root,subdir,files in sorted(os.walk(dirpath)) 
        for fn in sorted(files) 
        if fn.endswith(postfix)]
results = [p.get() for p in proc]

#%%
# geomcols=np.arange(8)
# wellcols=np.arange(12)
# CASSMGeom=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
#                         sheet_name='Geode',header=0,usecols=geomcols)
# CASSMsrc=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
#                         sheet_name='CASSM_Channels',header=0,usecols=geomcols)
# WellGeom=pd.read_excel('/home/spri902/Collab_metadata/Well_Points.xlsx',header=0,usecols=wellcols)
# fig=plt.figure()
# ax = fig.add_subplot(projection='3d')
# s1=ax.scatter(WellGeom.iloc[:,3:9:6]/3.28084,WellGeom.iloc[:,4:10:6]/3.28084,WellGeom.iloc[:,5:11:6]/3.28084,\
#             label='Well Points',marker='.',color='black',s=2)
# s2=ax.scatter(CASSMGeom.iloc[12:23,5],CASSMGeom.iloc[12:23,6],CASSMGeom.iloc[12:23,7],label='OT hydro',marker='o',color='blue')
# s3=ax.scatter(CASSMGeom.iloc[69:78:3,5],CASSMGeom.iloc[69:78:3,6],CASSMGeom.iloc[69:78:3,7],label='OT accels',marker='^',color='green')
# s4=ax.scatter(CASSMGeom.iloc[60:69:3,5],CASSMGeom.iloc[60:69:3,6],CASSMGeom.iloc[60:69:3,7],label='OB accels',marker='^',color='green')
# s5=ax.scatter(CASSMGeom.iloc[51:60:3,5],CASSMGeom.iloc[51:60:3,6],CASSMGeom.iloc[51:60:3,7],label='PST accels',marker='^',color='green')
# s6=ax.scatter(CASSMsrc.iloc[0:6,5],CASSMsrc.iloc[0:6,6],CASSMsrc.iloc[0:6,7],label='OB Sources',marker='*',color='red')
# s7=ax.scatter(CASSMsrc.iloc[6:12,5],CASSMsrc.iloc[6:12,6],CASSMsrc.iloc[6:12,7],label='PST Sources',marker='*',color='red')
# s8=ax.scatter(CASSMsrc.iloc[12:17,5],CASSMsrc.iloc[12:17,6],CASSMsrc.iloc[12:17,7],label='PSB Sources',marker='*',color='magenta')
# s9=ax.scatter(CASSMsrc.iloc[17,5],CASSMsrc.iloc[17,6],CASSMsrc.iloc[17,7],label='PDT Source',marker='*',color='red')
# s10=ax.scatter(CASSMGeom.iloc[:12,5],CASSMGeom.iloc[:12,6],CASSMGeom.iloc[:12,7],label='PDB hydro',marker='o',color='blue')
# s11=ax.scatter(CASSMGeom.iloc[24:30,5],CASSMGeom.iloc[24:30,6],CASSMGeom.iloc[24:30,7],label='PDT accels',marker='^',color='green')
# s12=ax.scatter(CASSMGeom.iloc[30:36,5],CASSMGeom.iloc[30:36,6],CASSMGeom.iloc[30:36,7],label='PDB accels',marker='^',color='green')
# s13=ax.scatter(CASSMGeom.iloc[36:39,5],CASSMGeom.iloc[36:39,6],CASSMGeom.iloc[36:39,7],label='PDT accels',marker='^',color='green')
# s14=ax.scatter(CASSMGeom.iloc[39:42,5],CASSMGeom.iloc[39:42,6],CASSMGeom.iloc[39:42,7],label='PDB accels',marker='^',color='green')
# s15=ax.scatter(CASSMGeom.iloc[42:51,5],CASSMGeom.iloc[42:51,6],CASSMGeom.iloc[42:51,7],label='PSB accels',marker='^',color='green')

# ax.set_xlabel('Easting')
# ax.set_ylabel('Northing')
# ax.set_zlabel('Elevation')
# ax.legend()
# plt.show()