#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 12:46:03 2022

simple script to do some CASSM data exploration and QC
look at single source single receiver data through time

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
from phasepickers.parker_cfs import CASSM_chanDecode, src_rec_dist
from datetime import datetime as datetime
from scipy.spatial.distance import cdist
import pickle

geomcols=np.arange(8)
wellcols=np.arange(12)
CASSMGeom=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='Geode',header=0,usecols=geomcols)
CASSMsrc=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='CASSM_Channels',header=0,usecols=geomcols)
CASSMsrc.columns = [c.replace(' ','_') for c in CASSMsrc.columns]
WellGeom=pd.read_excel('/home/spri902/Collab_metadata/Well_Points.xlsx',header=0,usecols=wellcols)
# meqcat = pd.read_csv('~/EGS_Collab/MEQ/catalog_manual_all_v190312.csv',header=0)
#%%

fig=plt.figure()
ax = fig.add_subplot(projection='3d')
s1=ax.scatter(WellGeom.iloc[1238:,3:9:6]/3.28084,WellGeom.iloc[1238:,4:10:6]/3.28084,WellGeom.iloc[1238:,5:11:6]/3.28084,\
            label='Monitoring Wells',marker='.',color='black',s=2)
s2=ax.scatter(CASSMGeom.iloc[12:23,5],CASSMGeom.iloc[12:23,6],CASSMGeom.iloc[12:23,7],label='OT hydro',marker='o',color='blue')
s3=ax.scatter(CASSMGeom.iloc[69:78:3,5],CASSMGeom.iloc[69:78:3,6],CASSMGeom.iloc[69:78:3,7],label='OT accels',marker='^',color='green')
s4=ax.scatter(CASSMGeom.iloc[60:69:3,5],CASSMGeom.iloc[60:69:3,6],CASSMGeom.iloc[60:69:3,7],label='OB accels',marker='^',color='green')
s5=ax.scatter(CASSMGeom.iloc[51:60:3,5],CASSMGeom.iloc[51:60:3,6],CASSMGeom.iloc[51:60:3,7],label='PST accels',marker='^',color='green')
s6=ax.scatter(CASSMsrc.iloc[0:6,5],CASSMsrc.iloc[0:6,6],CASSMsrc.iloc[0:6,7],label='OB Sources',marker='*',color='red')
s7=ax.scatter(CASSMsrc.iloc[6:12,5],CASSMsrc.iloc[6:12,6],CASSMsrc.iloc[6:12,7],label='PST Sources',marker='*',color='red')
s8=ax.scatter(CASSMsrc.iloc[12:17,5],CASSMsrc.iloc[12:17,6],CASSMsrc.iloc[12:17,7],label='PSB Sources',marker='*',color='magenta')
s9=ax.scatter(CASSMsrc.iloc[17,5],CASSMsrc.iloc[17,6],CASSMsrc.iloc[17,7],label='PDT Source',marker='*',color='red')
s10=ax.scatter(CASSMGeom.iloc[:12,5],CASSMGeom.iloc[:12,6],CASSMGeom.iloc[:12,7],label='PDB hydro',marker='o',color='blue')
s11=ax.scatter(CASSMGeom.iloc[24:30,5],CASSMGeom.iloc[24:30,6],CASSMGeom.iloc[24:30,7],label='PDT accels',marker='^',color='green')
s12=ax.scatter(CASSMGeom.iloc[30:36,5],CASSMGeom.iloc[30:36,6],CASSMGeom.iloc[30:36,7],label='PDB accels',marker='^',color='green')
s13=ax.scatter(CASSMGeom.iloc[36:39,5],CASSMGeom.iloc[36:39,6],CASSMGeom.iloc[36:39,7],label='PDT accels',marker='^',color='green')
s14=ax.scatter(CASSMGeom.iloc[39:42,5],CASSMGeom.iloc[39:42,6],CASSMGeom.iloc[39:42,7],label='PDB accels',marker='^',color='green')
s15=ax.scatter(CASSMGeom.iloc[42:51,5],CASSMGeom.iloc[42:51,6],CASSMGeom.iloc[42:51,7],label='PSB accels',marker='^',color='green')
s16=ax.scatter(WellGeom.iloc[:635,3:9:6]/3.28084,WellGeom.iloc[:635,4:10:6]/3.28084,WellGeom.iloc[:635,5:11:6]/3.28084,\
            label='Injector',marker='.',color='green',s=2)
s16=ax.scatter(WellGeom.iloc[635:1238,3:9:6]/3.28084,WellGeom.iloc[635:1238,4:10:6]/3.28084,WellGeom.iloc[635:1238,5:11:6]/3.28084,\
            label='Producer',marker='.',color='red',s=2)
# s17=ax.scatter(meqcat.x[:730],meqcat.y[:730],meqcat.z[:730],'+',color='orange')

ax.set_xlabel('Easting')
ax.set_ylabel('Northing')
ax.set_zlabel('Elevation')
# ax.legend()
plt.show()

#%%
# get file info and dates from directories
directory = sorted(os.walk('/home/spri902/mayCASSM'))
dates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]
xlims = mdates.date2num(dates)

# src=[0] # OB source shallowest
# # OT hydrophone and accels listed by channel number ** not index **
# recVec = np.array((13, 14, 15, 16, 17, 18, 19, 20, 21,
#           22, 23, 70, 71, 72, 73, 74, 75, 76, 77, 78))
# # channels for OT hydrophones and accels
# hydChan=slice(12,23)
# accChan=slice(69,78)

# src=[5] # OB source deepest
# # PST and PDT accels listed by channel number ** not index **
# recVec = np.array((25,26,27,28,29,30,37,38,39,52,53,54,55,56,57,58,59,60))
# # channels for PDT and PST hydrophones and accels

# accChan1=slice(24,30) # index not channel numbers
# accChan2=slice(36,39) # index not channel numbers
# accChan3=slice(51,60) # index not channel numbers

# PST src and channels for PDB hydrophones and accels
src=[0] # PST source
# hydChan=slice(0,12) # PDB hydrophones
# accChan=slice(30,36)# PDB accels
# ac2Chan=slice(39,42)# PDB accels
accChan = slice(54,57)
# PDB hydrophone and accels listed by channel number ** not index **
# recVec = np.array((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
#                     31, 32, 33, 34, 35, 36, 40, 41, 42))
recVec = np.array((55,56,57))
# pull src and rec locations
srcloc = CASSMsrc.loc[src[0],'Easting_m':'Elev_m'].values
srcloc = srcloc.reshape(-1,3).astype(float)
recloc = CASSMGeom.loc[recVec-1,'Easting(m)':'Elev(m)'].values
recInd = np.arange(len(recVec))
# only need to grab part of the signal (don't need all 9600 samples)
siglen = 2000
samples = np.arange(siglen)

OThydata = np.empty((len(dates),siglen,len(recVec)))
cnt=0
low=3000
high=8000

#%% Build the single src single rec gather through time and save as numpy array 
for path,subdir,files in sorted(os.walk('/home/spri902/mayCASSM')):
    for filename in sorted(files):
        if filename.endswith('.dat'):
            sthead = _read_seg2(os.path.join(path,filename), headonly=True)
            times = sthead[0].times()
            srcCh = CASSM_chanDecode(t=times, enc=sthead[93].data)
            if srcCh in src:
                srcID = CASSMsrc[CASSMsrc["CytekCh"] == srcCh]
                dt = sthead[0].stats.delta
                fs = np.floor(sthead[0].stats.sampling_rate)
                st = sthead[0].stats.starttime
                et = sthead[0].stats.endtime
                sthead.trim(st, st + (dt*1999))
                # OThy = sthead[hydChan]+sthead[accChan]+sthead[ac2Chan]
                OThy = sthead[accChan]
                                                              
                OThy = OThy.filter('bandpass', freqmin=low,freqmax=high,
                                    corners=4, zerophase=True)
                OThydata[cnt,:,:] = np.stack(tr.data for tr in OThy).T
                
                cnt += 1
                
            else:
                continue
os.chdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/src0_PSTacc')            
with open('OB0data.npy','wb') as f:
    np.save(f,OThydata)


#%% make plots of single src single rec gather through time with pcolormesh
# load data if created and plot
os.chdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/src0_PSTacc')
OThydata = np.load('OB0data.npy')
for i in recInd:    
    fig1,ax = plt.subplots()
    recID = CASSMGeom[CASSMGeom["Channel"] == recVec[i]].Sensor.item()
    dist = np.round(src_rec_dist(srcloc,recloc),2)
    vm = np.percentile(OThydata[:,:,i],99)
    #img = ax.imshow(OThydata[:,:,i].T,extent=[xlims[0],xlims[-1],1999,0],aspect='auto',cmap='RdBu',vmin=-vm,vmax=vm)
    img = ax.pcolormesh(xlims,samples,OThydata[:,:,i].T,cmap='RdBu')
    plt.gca().invert_yaxis()
    ax.xaxis_date()
    plt.colorbar(img)
    plt.gcf().autofmt_xdate()
    plt.title(f'CASSM Source: {src[0]} Receiver: {recID} S_R_dist: {dist[0][i]}')
    plt.axvline(dates[42],color='green')
    plt.axvline(dates[50],color='red')
    plt.axvline(dates[96],color='green')
    plt.axvline(dates[106],color='red')
    plt.axvline(dates[224],color='green')
    plt.axvline(dates[232],color='red')
    plt.axvline(dates[352],color='green')
    plt.axvline(dates[368],color='red')
    plt.axvline(dates[472],color='green')
    plt.axvline(dates[475],color='red')
    plt.axvline(dates[507],color='green')
    plt.axvline(dates[512],color='red')
    pickle.dump(ax, open(f'CASSM_Source_{src[0]}_Receiver_{recID}','wb'))
    plt.close()
    
for i in recInd:
    recID = CASSMGeom[CASSMGeom["Channel"] == recVec[i]].Sensor.item()
    fig = pickle.load(open(f'CASSM_Source_{src[0]}_Receiver_{recID}','rb'))


#%% Try autocorrelations
# def autocorr(x):
#     result = np.correlate(x, x, mode='full')
#     return result[result.size//2:]

# ac = autocorr(OThydata[597,:,16])


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    