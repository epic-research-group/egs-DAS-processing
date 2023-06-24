#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 12:44:41 2022

@author: spri902
"""

import matplotlib.pyplot as plt
%matplotlib qt5
import os
import numpy as np
import pandas as pd
import CASSM_decoder as decoder
from obspy.io.seg2.seg2 import _read_seg2
from scipy.signal import *
import dtwFuncs as dtw



# read in the well geometry and CASSM src / rec geometry
geomcols = np.arange(8)
wellcols = np.arange(12)

CASSMrec = pd.read_excel(
    '~/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='Geode',header=0,usecols=geomcols)
CASSMsrc = pd.read_excel(
    '~/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='CASSM_Channels',header=0,usecols=geomcols)
WellGeom = pd.read_excel(
    '/home/spri902/Collab_metadata/Well_Points.xlsx',header=0,usecols=wellcols)

src = [1,2] # OB sources
rec = [22,23,24]
recChan = slice(22,25) # OT receivers
ptd = '/home/spri902/EGS_Collab/maystim'

numDates = len(sorted(os.walk(ptd))[0][1]) # num of dates to pull data from
#setup the list of lists for the data
hydata = [[[] for i in range(numDates)] for i in range(len(src))]
t = [[[] for i in range(numDates)] for i in range(len(src))]
id1 = 0
id2 = 0
#%% put data into list of lists rows are numSrcs cols are time directories

# data sorted from early to late (ie column 1 is early last column is late)

for path,subdir,files in sorted(os.walk(ptd)):
    for filename in sorted(files):
        if filename.endswith(".dat"):
            sthead = _read_seg2(os.path.join(path,filename), headonly=True)
            times = sthead[0].times()
            srcCh = decoder.CASSM_chanDecode(t=times, enc=sthead[93].data)
            srcID = CASSMsrc[CASSMsrc["CytekCh"] == srcCh]
            dt = sthead[0].stats.delta
            fs = np.floor(sthead[0].stats.sampling_rate)
            nyq = np.floor(0.5*(1/dt))
            if srcCh in src:
                if srcCh == src[0]:
                    st = sthead[0].stats.starttime
                    et = sthead[0].stats.endtime
                    sthead.trim(st, st + (dt*1439))
                    hydata[0][id1] = np.stack(tr.data for tr in sthead[recChan]).T
                    t[0][id1]= np.stack(tr.times() for tr in sthead[recChan]).T
                    # print(sthead[recChan].stats.seg2['RAW_RECORD'])
                    id1 += 1
                else:
                    st = sthead[0].stats.starttime
                    et = sthead[0].stats.endtime
                    sthead.trim(st, st + (dt*1439))
                    hydata[1][id2] = np.stack(tr.data for tr in sthead[recChan]).T
                    t[1][id2]= np.stack(tr.times() for tr in sthead[recChan]).T
                    # print(sthead[recChan].stats.seg2['RAW_RECORD'])
                    id2 += 1
            else: 
                continue
#%% save the data as a pickle file
# dtw.save_file(hydata,'src_1_2_rec_22')
#hydata = dtw.load_object('src_1_2_rec_22.pickle')
#%%
npts = len(hydata[0][0])
trc_data = np.zeros((npts,numDates))
tim_data = np.zeros((npts,numDates))
taper_win = tukey(npts,alpha=0.05)

Wn = [2000,3000] / np.array(nyq) # [Hz] low and high freqs for the filter
sos = butter(4,Wn,btype='bandpass',output='sos')

n_src = len(src)
n_rec = len(rec)

maxLag=200
b=4

t_start = 0.01
n_start = np.argmin(abs(tim_data[:,1] - t_start))
#dv_v = [[[ [] for k in range(n_rec)] for j in range(n_src)] for i in range(numDates)]
dv_v = np.zeros((n_src,n_rec,numDates))

for isrc in range(n_src):
    for irec in range(n_rec):
        for iday in range(numDates):
            recData = hydata[isrc][iday]
            if recData.size > 0:
                trc_data[:,iday] = detrend(recData[:,irec]*taper_win)
                tim_data[:,iday] = t[isrc][iday][:,0]
        
        trc_data = sosfiltfilt(sos,trc_data)
        
        for iday in range(numDates):
            err = dtw.computeErrFunc(trc_data[:,iday], trc_data[:,0], npts,lag=maxLag)
            
            #accumulate error in forward direction
            direction = 1 # direction to accumulate error (1=forward -1=backwrd)
            # now do forward accumulation to make distance fcn
            dist = dtw.accumulateErrFunc(direction, err, npts, maxLag, b)
            # find shifts in backward direction
            stbar = dtw.backtrackDistFunc(-direction, dist, err, -maxLag, b)
            stbarTime = stbar * dt # convert from samples to time
            # linear polynomial fitting
            p = np.polyfit(tim_data[n_start:,0], stbar[n_start:],1)
            dv_v = [isrc,irec,iday] = -p[0]
            
            
            
            
            
            
            
            
            
            
            
            
#%%
fig=plt.figure()
ax = fig.add_subplot(projection='3d')
s1=ax.scatter(WellGeom.iloc[:,3:9:6]/3.28084,WellGeom.iloc[:,4:10:6]/3.28084,WellGeom.iloc[:,5:11:6]/3.28084,\
           label='Well Points',marker='.',color='black',s=2)
s2=ax.scatter(CASSMrec.iloc[12:23,5],CASSMrec.iloc[12:23,6],CASSMrec.iloc[12:23,7],label='OT hydro',marker='o',color='blue')
s3=ax.scatter(CASSMrec.iloc[69:78:3,5],CASSMrec.iloc[69:78:3,6],CASSMrec.iloc[69:78:3,7],label='OT accels',marker='^',color='green')
s4=ax.scatter(CASSMrec.iloc[60:69:3,5],CASSMrec.iloc[60:69:3,6],CASSMrec.iloc[60:69:3,7],label='OB accels',marker='^',color='green')
s5=ax.scatter(CASSMrec.iloc[51:60:3,5],CASSMrec.iloc[51:60:3,6],CASSMrec.iloc[51:60:3,7],label='PST accels',marker='^',color='green')
s6=ax.scatter(CASSMsrc.iloc[0:6,5],CASSMsrc.iloc[0:6,6],CASSMsrc.iloc[0:6,7],label='OB Sources',marker='*',color='red')
s7=ax.scatter(CASSMsrc.iloc[6:12,5],CASSMsrc.iloc[6:12,6],CASSMsrc.iloc[6:12,7],label='PST Sources',marker='*',color='red')
s8=ax.scatter(CASSMsrc.iloc[12:17,5],CASSMsrc.iloc[12:17,6],CASSMsrc.iloc[12:17,7],label='PSB Sources',marker='*',color='magenta')
s9=ax.scatter(CASSMsrc.iloc[17,5],CASSMsrc.iloc[17,6],CASSMsrc.iloc[17,7],label='PDT Source',marker='*',color='red')
s10=ax.scatter(CASSMrec.iloc[:12,5],CASSMrec.iloc[:12,6],CASSMrec.iloc[:12,7],label='PDB hydro',marker='o',color='blue')

ax.set_xlabel('Easting')
ax.set_ylabel('Northing')
ax.set_zlabel('Elevation')
ax.legend()
plt.show()