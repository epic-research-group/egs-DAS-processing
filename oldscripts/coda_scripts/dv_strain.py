#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:02:06 2022

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
# %% load in accelerometer Xcorr coeffs and strain rate data
cassmdir = 'src9_PDBhyds_accs'
os.chdir(f'/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/{cassmdir}')

stretchdata = np.load('PST9_stretch_16.npz')
dtArray = stretchdata['arr_0']
tSamp = stretchdata['arr_1']
cArray = stretchdata['arr_2']

ccorrdata = np.load('PST9_ccorr.npy')
ccorr = ccorrdata['arr_0']

os.chdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/wellPDB')
datFiles = sorted(os.listdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/wellPDB'))
df_full = pd.concat((pd.read_csv(f,header=0,sep=',') for f in datFiles if f.endswith('allchan.txt')),axis=0)

directory = sorted(os.walk('/home/spri902/mayCASSM'))
CASSMdates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]

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
    
dasdnums = mdates.date2num(DASdates)
cassmdnums = mdates.date2num(CASSMdates)

# df_full = detrend(df_full,axis=0,type='constant')
DASint = resample(df_full,len(cassmdnums[:-2]))
dasd = pd.date_range(start=DASdates[0],end=DASdates[-1],periods=len(cassmdnums[:-2]))
df_full = df_full*116.0
chans=np.linspace(0,df_full.shape[1] - 1,df_full.shape[1]).astype(int)
stimbeg = [96, 224, 352, 472, 507]
stimfin = [106, 232, 368, 475, 512]
stimbegLines = itemgetter(*stimbeg)(CASSMdates)
stimfinLines = itemgetter(*stimfin)(CASSMdates)

# %%
os.chdir('/home/spri902/EGS_Collab/4850/stimflow/')
injFiles = sorted(os.listdir('/home/spri902/EGS_Collab/4850/stimflow/'))
injDat = pd.concat((pd.read_csv(f,header=1,usecols=[0,24],\
                           parse_dates = [0],infer_datetime_format=True) \
                    for f in injFiles if f.endswith('.csv')),axis=0)
injDat.rename(columns={'hh:mm:ss':'date','psig.9':'psig'},inplace=True)
# tmp = pd.read_csv(injFiles[3],header=1,usecols=[0,16],\
                          # parse_dates=[0],infer_datetime_format=True)
# tmp.rename(columns={'hh:mm:ss':'date','psig.2':'psig'},inplace=True)
# injDat = pd.concat((injDat,tmp),axis=0)
injDat.reset_index(drop = True,inplace = True)
injDat.set_index('date',inplace=True)
# %%
fig, (ax,ax1) = plt.subplots(2,1,sharex=True,figsize=(12,12))
fig.subplots_adjust(right=0.75)

img1 = ax.pcolormesh(cassmdnums[:-2],tSamp,dtArray[:-1,:].T*100,cmap='RdBu')
ax.invert_yaxis()
ax.set_ylim(0.04,0.005)

ax1 = injDat.loc[:,("psig")].plot()
ax1.set_ylabel('Injection Pressure (psig)')

ax2=ax.twinx()
ax3=ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
ax.xaxis_date()
ax3.spines.right.set_position(("axes",1.1))
ax2.set_ylim(-1,0)
ax3.set_ylim(0,1)
ax2.set_ylabel('DAS strain rate')
ax3.set_ylabel('CASSM Xcorr')

img2 = ax2.plot(dasd,DASint[:,30],"b",label="DAS Strain Rate")
img3 = ax3.plot(CASSMdates[1:],ccorr[16,:],"orange",label='Zero-Lag Corrcoeff for PDB ')
ax2.yaxis.label.set_color(img2[0].get_color())
ax3.yaxis.label.set_color(img3[0].get_color())
cbar = plt.colorbar(img1,orientation="horizontal")
cbar.set_label('dV in % velocity change')
tkw = dict(size=4, width=1.5)
ax2.tick_params(axis='y', colors=img2[0].get_color(), **tkw)
ax3.tick_params(axis='y', colors=img3[0].get_color(), **tkw)
ax1.tick_params(axis='x', **tkw)

l1 = ax.axvspan(stimbegLines[1],stimfinLines[1],alpha=0.5,color='grey',label='Injection')
l2 = ax.axvspan(stimbegLines[2],stimfinLines[2],alpha=0.5,color='grey',label='Injection')
l3 = ax.axvspan(stimbegLines[3],stimfinLines[3],alpha=0.5,color='grey',label='Injection')
l4 = ax.axvspan(stimbegLines[4],stimfinLines[4],alpha=0.5,color='grey',label='Injection')
lns = img2+img3+[l1]

labls = [l.get_label() for l in lns]
ax.legend(lns,labls, loc='lower right')
plt.title('DAS strain rate and Acclerometer Xcorr plotted on top of Trace Stretching image')

# ax1[1].annotate(
#     'Injection 1',
#     xy=(0, 1), xycoords='data',
#     xytext=(-50, 30), textcoords='offset points',
#     arrowprops=dict(arrowstyle="->"))

plt.show()
# %%
fig,ax = plt.subplots()
# sc = ax.scatter(CASSMint,PDBdat[29,:],c=dasdnums,marker='.',cmap='PiYG')
sc = ax.scatter(ccorr[16,:-1],DASint[:,30],c=dasd,marker='.',cmap='PiYG')
# sc = ax.scatter(ccorr[16,:-1],DASint,c=cassmdnums[:-2],marker='.',cmap='PiYG')
loc = mdates.AutoDateLocator()
fig.colorbar(sc,ticks=loc,format=mdates.AutoDateFormatter(loc))
plt.xlabel('CorrCoeff values')
plt.ylabel('DAS Strain Rate')
plt.title('DAS strain rate vs CorrCoeff values')
plt.xlim([1,0.1])
# plt.xlim([0.04,0])