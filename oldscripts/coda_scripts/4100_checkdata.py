#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:45:02 2022

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
#%%


Geom=pd.read_csv('/home/spri902/Collab_metadata/4100ChannelMappingGeode.csv',\
                      header=[0])

CASSMsrc=pd.read_csv('/home/spri902/Collab_metadata/4100ChannelMappingSource.csv',\
                        header=[0])

#%%
# get file info and dates from directories
directory = sorted(os.walk('/shared/Collab4100/may2022stim'))
dates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]
xlims = mdates.date2num(dates)


# AML sources and hydrophones in TS
src=[8] # AML source
hydChan=slice(48,60) # TS hydrophones
srcID = CASSMsrc['Borehole'][src[0]]
# TS hydrophones listed by channel number ** not index **
recVec = np.arange(49,61) # remember this is by Geode num so Index + 1
# pull src and rec locations
# srcloc = CASSMsrc.loc[src[0],'Easting_m':'Elev_m'].values
# srcloc = srcloc.reshape(-1,3).astype(float)
# recloc = CASSMGeom.loc[recVec-1,'Easting(m)':'Elev(m)'].values
recInd = np.arange(len(recVec))
# only need to grab part of the signal (don't need all 9600 samples)
siglen = 3000
samples = np.arange(siglen)

TShydata = np.empty((len(dates),siglen,len(recVec)))
cnt=0
low=1000
# high=10000
#%% Build the single src single rec gather through time and save as numpy array 
for path,subdir,files in sorted(os.walk('/shared/Collab4100')):
    for filename in sorted(files):
        if filename.endswith('.dat'):
            sthead = _read_seg2(os.path.join(path,filename), headonly=True)
            times = sthead[0].times()
            srcCh = CASSM_chanDecode(t=times, enc=sthead[73].data)
            if srcCh in src:
                
                TShy = sthead[hydChan].copy()
                TShy.detrend('constant') # mean of data is subtracted
                TShy.detrend('simple') # subtract linear fcn defined by 1st/last samples of the trace     
                dt = TShy[0].stats.delta
                fs = np.floor(TShy[0].stats.sampling_rate)
                st = TShy[0].stats.starttime
                et = TShy[0].stats.endtime
                TShy.trim(st, st + (dt*2999))
                
                                                              
                TShy = TShy.filter('highpass', freq=low,corners=4, zerophase=True)
                TShydata[cnt,:,:] = np.stack(tr.data for tr in TShy).T
                
                cnt += 1
                
            else:
                continue
            
with open('/home/spri902/EGS_Collab/4100/results/marchstim/TS10data.npy','wb') as f:
    np.save(f,TShydata)
    #%% make plots of single src single rec gather through time with pcolormesh
    # load data if created and plot
   TShydata = np.load('TS8data.npy')
    for i in recInd:    
        fig1,ax = plt.subplots()
        recID = Geom[Geom["Geode"] == recVec[i]].Component.item()
        # dist = np.round(src_rec_dist(srcloc,recloc),2)
        vm = np.percentile(TShydata[:,:,i],99)
        #img = ax.imshow(OThydata[:,:,i].T,extent=[xlims[0],xlims[-1],1999,0],aspect='auto',cmap='RdBu',vmin=-vm,vmax=vm)
        img = ax.pcolormesh(xlims,samples,TShydata[:,:,i].T,cmap='RdBu',vmin=-vm,vmax=vm)
        plt.gca().invert_yaxis()
        ax.xaxis_date()
        plt.colorbar(img)
        plt.gcf().autofmt_xdate()
        plt.title(f'CASSM Source: {srcID} {src[0]} Receiver: {recID} ')
        pickle.dump(ax, open(f'CASSM_Source_{src[0]}_Receiver_{recID}','wb'))
        plt.close()
        
    for i in recInd:
        recID = Geom[Geom["Geode"] == recVec[i]].Component.item()
        fig = pickle.load(open(f'CASSM_Source_{src[0]}_Receiver_{recID}','rb'))