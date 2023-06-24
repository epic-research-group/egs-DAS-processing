#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 12:44:41 2022

@author: spri902
"""

import matplotlib
%matplotlib
import matplotlib.pyplot as plt 
import os
import numpy as np
import pandas as pd
from datetime import datetime as datetime
from phasepickers.parker_cfs import CASSM_chanDecode, src_rec_dist
from coda_scripts.dtwFuncs import *
# import cupy as cp

#%%


# read in the well geometry and CASSM src / rec geometry
geomcols = np.arange(8)
wellcols = np.arange(12)

CASSMrec = pd.read_excel(
    '~/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='Geode',header=0,usecols=geomcols)
CASSMsrc = pd.read_excel(
    '~/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='CASSM_Channels',header=0,usecols=geomcols)
CASSMGeom=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                            sheet_name='Geode',header=0,usecols=geomcols)
WellGeom = pd.read_excel(
    '/home/spri902/Collab_metadata/Well_Points.xlsx',header=0,usecols=wellcols)

# PST src and channels for PDB hydrophones and accels
src=[11] # PST source
hydChan=slice(0,12) # PDB hydrophones
accChan=slice(30,36)# PDB accels
ac2Chan=slice(39,42)# PDB accels
# PDB hydrophone and accels listed by channel number ** not index **
recVec = np.array((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                    31, 32, 33, 34, 35, 36, 40, 41, 42))
# pull src and rec locations
srcloc = CASSMsrc.loc[src[0],'Easting_m':'Elev_m'].values
srcloc = srcloc.reshape(-1,3).astype(float)
recloc = CASSMGeom.loc[recVec-1,'Easting(m)':'Elev(m)'].values
recInd = np.arange(len(recVec))
# only need to grab part of the signal (don't need all 9600 samples)
siglen = 2000
samples = np.arange(siglen)
dat_path = '/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/src9_PDBhyds_accs/'
fn = 'PST9data.npy'

seisdata = np.load(os.path.join(dat_path,fn))
# src = [1,2] # OB sources
# rec = [22,23,24]
# recChan = slice(22,25) # OT receivers
# ptd = '/home/spri902/EGS_Collab/maystim'

# numDates = len(sorted(os.walk(ptd))[0][1]) # num of dates to pull data from
# #setup the list of lists for the data
# hydata = [[[] for i in range(numDates)] for i in range(len(src))]
# t = [[[] for i in range(numDates)] for i in range(len(src))]
# id1 = 0
# id2 = 0
#%% put data into list of lists rows are numSrcs cols are time directories

# data sorted from early to late (ie column 1 is early last column is late)

# for path,subdir,files in sorted(os.walk(ptd)):
#     for filename in sorted(files):
#         if filename.endswith(".dat"):
#             sthead = _read_seg2(os.path.join(path,filename), headonly=True)
#             times = sthead[0].times()
#             srcCh = decoder.CASSM_chanDecode(t=times, enc=sthead[93].data)
#             srcID = CASSMsrc[CASSMsrc["CytekCh"] == srcCh]
#             dt = sthead[0].stats.delta
#             fs = np.floor(sthead[0].stats.sampling_rate)
#             nyq = np.floor(0.5*(1/dt))
#             if srcCh in src:
#                 if srcCh == src[0]:
#                     st = sthead[0].stats.starttime
#                     et = sthead[0].stats.endtime
#                     sthead.trim(st, st + (dt*1439))
#                     hydata[0][id1] = np.stack(tr.data for tr in sthead[recChan]).T
#                     t[0][id1]= np.stack(tr.times() for tr in sthead[recChan]).T
#                     # print(sthead[recChan].stats.seg2['RAW_RECORD'])
#                     id1 += 1
#                 else:
#                     st = sthead[0].stats.starttime
#                     et = sthead[0].stats.endtime
#                     sthead.trim(st, st + (dt*1439))
#                     hydata[1][id2] = np.stack(tr.data for tr in sthead[recChan]).T
#                     t[1][id2]= np.stack(tr.times() for tr in sthead[recChan]).T
#                     # print(sthead[recChan].stats.seg2['RAW_RECORD'])
#                     id2 += 1
#             else: 
#                 continue
#%% save the data as a pickle file
# dtw.save_file(hydata,'src_1_2_rec_22')
#hydata = dtw.load_object('src_1_2_rec_22.pickle')
#%%
directory = sorted(os.walk('/home/spri902/mayCASSM'))
CASSMdates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]

npts = seisdata.shape[1]
maxLag=30
b=7
dt=1/48000
epsilon=[]
epsTime=[]
t_start = 0.004
CASSMtimes = np.r_[0:seisdata.shape[1]]*dt
# n_start = np.argmin(abs(tim_data[:,1] - t_start))
#dv_v = [[[ [] for k in range(n_rec)] for j in range(n_src)] for i in range(numDates)]
# dv_v = np.zeros((n_src,n_rec,numDates))
epsilon = np.zeros((siglen,0))
epsTime = np.zeros((siglen,0))
for i in range(1,seisdata.shape[0]):
    u0 = seisdata[0,:,16]
    u1 = seisdata[i,:,16]

    err = computeErrFunc(u0,u1, npts,lag=maxLag)
    
    #accumulate error in forward direction
    direction = 1 # direction to accumulate error (1=forward -1=backwrd)
    # now do forward accumulation to make distance fcn
    dist = accumulateErrFunc(direction, err, npts, maxLag, b)
    # find shifts in backward direction
    stbar = backtrackDistFunc(-direction, dist, err, -maxLag, b).reshape((siglen,1))
    stbarTime = stbar * dt # convert from samples to time
    epsilon = np.append(epsilon,stbar,axis=1)
    epsTime = np.append(epsTime,stbarTime,axis=1)
    print(i)
    # linear polynomial fitting
    # p = np.polyfit(tim_data[n_start:,0], stbar[n_start:],1)
    # dv_v = [isrc,irec,iday] = -p[0]
    
#%%
chans=np.linspace(0,len(epsilon) ,len(epsilon)).astype(int)            
fig1,ax = plt.subplots()

vm = np.nanpercentile(epsTime,99)
# if vm < 0:
img=ax.pcolormesh(CASSMdates[1:],CASSMtimes,epsTime,cmap='RdBu',vmin=-vm,vmax=vm)
# if vm > 0:
    # img=ax.pcolormesh(dates,chans,df_all,cmap='RdBu',vmin=vm/2,vmax=vm*2)
plt.gca().invert_yaxis()
#ax.hlines(min(xlims),max(xlims),color='k')
#ax.text(xlims[50],400,"OT well", color='black', ha="center", va="center", bbox =dict(facecolor='none', boxstyle='square, pad=0.3', lw=2))
ax.xaxis_date()

#ax2.set_ylabels(wn)
plt.colorbar(img,orientation="horizontal")
#plt.gcf().autofmt_xdate()
plt.title(f'iDAS data {CASSMdates[0]} to {CASSMdates[-1]}')            
            
            
    