#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:11:05 2022

@author: spri902
"""

import os
import math
import pandas as pd
import numpy as np
import time
import matplotlib.dates as mdates
from datetime import datetime as datetime
import matplotlib
import matplotlib.pyplot as plt 
import pickle
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import pchip_interpolate
from scipy.signal import resample, detrend
from operator import itemgetter

# %% Read in DAS data

# and to read the pickle
# df_all = pd.read_pickle('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/wellPDB/maystim22_26_combined_n')
# df_r = pd.read_pickle('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/wellPDB/maystim22_26_combined_r')
df_full = pd.read_pickle('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/wellPDB/maystim22_26_combined_full')
dasScaler = 116.0
df_full = df_full.multiply(dasScaler)
df_full = detrend(df_full,axis=0,type='constant')
#mask the OT well data and remove column means from the remaining data
# OTchans  = np.r_[337:456]
# OBchans  = np.r_[518:635]
# PSTchans = np.r_[796:879]
# PSBchans = np.r_[923:1042]
# PDBchans = np.r_[1091:1210]
# PDTchans = np.r_[1271:1388]
 
# df_all = df_all - df_all.mean()

# df_all.mask((df_all.iloc[0:OTchans[-1]+1,:].notna()),inplace=True)
# df_all.mask((df_all.iloc[0:OTchans[0],:].notna()),inplace=True) 
# df_all.mask((df_all.iloc[OTchans[-1]+1:OBchans[0],:].notna()),inplace=True)  
# df_all.mask((df_all.iloc[OBchans[-1]+1:PSTchans[0],:].notna()),inplace=True) 
# df_all.mask((df_all.iloc[PSTchans[-1]+1:PSBchans[0],:].notna()),inplace=True)  
# df_all.mask((df_all.iloc[PSBchans[-1]+1:PDBchans[0],:].notna()),inplace=True) 
# df_all.mask((df_all.iloc[PDBchans[-1]+1:PDTchans[0],:].notna()),inplace=True)  
# df_all.mask((df_all.iloc[PDTchans[-1]+1:1728].notna()),inplace=True)

wn     = ['OT','OB','PST','PSB','PDB','PDT']
nfile_list = sorted(os.walk('/data1/parker/EGS_iDAS'))
nfile_list = nfile_list[1:]
#file_list = file_list[1:]
nfile_list = [group[2] for group in nfile_list]
nfile_list = [item for sublist in nfile_list for item in sublist]
# [file_list.append(f) for f in nfile_list]
fd = [name.split("_") for name in nfile_list]
fl = [fd[file][2].split(".") for file in range(len(fd))]
fl = [el[0] for el in fl]
DASdates = [datetime.strptime(d,'%y%m%d%H%M%S') for d in sorted(fl)]
ind2rem = [0, 90,   91, 257, 258, 1571, 1572, 3082, 3083, 5085, 5086, 5599, 5600, 5961, 5962, 7623, 7624, 8841, 8842, 9562]
# remove in reverse so that the indices remain in the correct order for removal
for index in sorted(ind2rem,reverse=True):
    del DASdates[index]
with open('DASdates.pkl','wb') as f:
    pickle.dump(DASdates,f)
# chans=np.linspace(0,1727,1728).astype(int)
# ch_bot = [396, 576, 837, 982, 1150, 1329]

# PDBdown = PDBchans[0:math.trunc(len(PDBchans)/2)]
# PDBup   = np.flipud(PDBchans[math.trunc(len(PDBchans)/2)+1:])
# PDBdat  = df_all.iloc[PDBdown,:].values + df_all.iloc[PDBup,:].values
# %% Read in CASSM data
directory = sorted(os.walk('/home/spri902/mayCASSM'))
CASSMdates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]
filepath = '/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/src9_PDBhyds_accs/'
datfile = 'PST9data.npy'
seisdata = np.load(os.path.join(filepath, datfile))
os.chdir(filepath)

stimbeg = [96, 224, 352, 472, 507]
stimfin = [106, 232, 368, 475, 512]
stimbegLines = itemgetter(*stimbeg)(CASSMdates)
stimfinLines = itemgetter(*stimfin)(CASSMdates)
stimvec=[]
[stimvec.append(np.r_[stimbeg[stim]:stimfin[stim]+1]) for stim in range(len(stimbeg))]
stimvec=np.concatenate(stimvec).ravel().tolist()
ccorr = np.load('PST9_ccorr.npy')['arr_0']


# %% Make new time axis
dasdnums = mdates.date2num(DASdates)
cassmdnums = mdates.date2num(CASSMdates)
CASSMint = pchip_interpolate(cassmdnums[:-2],ccorr[16,:-1],dasdnums)

DASint = resample(df_full,len(cassmdnums[:-2]))
dasd = pd.date_range(start=DASdates[0],end=DASdates[-1],periods=len(cassmdnums[:-2]))
chans=np.linspace(0,df_full.shape[1] - 1,df_full.shape[1]).astype(int)
# %%
fig,ax1 = plt.subplots()
ax2 = ax1.twinx()
# p1=ax1.plot(DASdates,PDBdat[29,:],label='Well PDB DAS strain rate')
p1=ax1.plot(dasd,DASint[:,30],label='Well PDB DAS strain rate')
# p1=ax1.plot(DASdates,df_all,label='Well PDB DAS strain rate')
p2=ax2.plot(CASSMdates[1:],ccorr[16,:],label='Zero-Lag Corrcoeff for PDB ',color='orange')
plt.vlines(stimbegLines,0,1,color='green')
plt.vlines(stimfinLines,0,1,color='red')
ax1.set_xlabel('Experiment Date')
ax1.set_ylabel('DAS Strain Rate')
ax2.set_ylabel('CASSM Receiver CorrCoeff values')
plt.title('PDB DAS strain rate CASSM receiver comparison')
lns = p1+p2
labls = [l.get_label() for l in lns]
ax1.legend(lns,labls, loc='lower right')
ax2.set_ylim([0.1,1])



# Upsampled CASSM corrcoeffs files to DAS files
fig,ax = plt.subplots()
# sc = ax.scatter(CASSMint,PDBdat[29,:],c=dasdnums,marker='.',cmap='PiYG')
sc = ax.scatter(ccorr[16,:-1],DASint,c=dasd,marker='.',cmap='PiYG')
# sc = ax.scatter(ccorr[16,:-1],DASint,c=cassmdnums[:-2],marker='.',cmap='PiYG')
loc = mdates.AutoDateLocator()
fig.colorbar(sc,ticks=loc,format=mdates.AutoDateFormatter(loc))
plt.xlabel('CorrCoeff values')
plt.ylabel('DAS Strain Rate')
plt.title('DAS strain rate vs CorrCoeff values')
plt.xlim([1,0.1])
# Downsampled DAS data to CASSM data
# fig,ax = plt.subplots()
# sc = ax.scatter(ccorr[16,:],DASint,c=cassmdnums[1:],marker='.',cmap='turbo')
# loc = mdates.AutoDateLocator()
# fig.colorbar(sc,ticks=loc,format=mdates.AutoDateFormatter(loc))
# plt.xlabel('CorrCoeff values')
# plt.ylabel('DAS Strain Rate')
# plt.title('DAS strain rate vs CorrCoeff values')

# %% plot image of entire well of DAS data


fig, ax = plt.subplots()
vm=np.nanpercentile(df_full,99)
if vm < 0:
    img=ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=2*vm,vmax=vm+abs(vm))
if vm > 0:
    img=ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=vm/2,vmax=vm*2)
# img = ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=-0.05,vmax=0.05)
plt.gca().invert_yaxis()
ax.xaxis_date()
plt.colorbar(img,orientation="horizontal")

fig, ax = plt.subplots()
vm=np.nanpercentile(df_full,99)
# if vm < 0:
#     img=ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=2*vm,vmax=vm+abs(vm))
# if vm > 0:
#     img=ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=vm/2,vmax=vm*2)
img = ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=-vm,vmax=vm)
plt.gca().invert_yaxis()
ax.xaxis_date()
plt.colorbar(img,orientation="horizontal")

df_strain = cumulative_trapezoid(DASint,axis=0,initial=0)
# df_strain = detrend(df_strain,axis = 0,type = 'linear')
# df_strain = detrend(df_strain,axis = 0,type = 'constant')

fig, ax = plt.subplots()
vm=np.nanpercentile(df_strain,99)
# if vm < 0:
#     img=ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=2*vm,vmax=vm+abs(vm))
# if vm > 0:
#     img=ax.pcolormesh(DASdates,chans,df_full.T,cmap='RdBu',vmin=vm/2,vmax=vm*2)
img = ax.pcolormesh(dasd,chans,df_strain.T,cmap='RdBu',vmin=-vm,vmax=vm)
plt.gca().invert_yaxis()
ax.xaxis_date()
plt.colorbar(img,orientation="horizontal")