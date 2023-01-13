#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 22:18:37 2021

@author: spri902
"""

#%%
import pickle
import numpy as np
import pandas as pd
%matplotlib
import matplotlib.pyplot as plt
import os
from obspy.io.seg2.seg2 import _read_seg2
import scaleogram as scg
from phasepickers import *
from scipy.signal import *
#from scipy import stats
#from scipy.optimize import curve_fit
#from tsmoothie.smoother import *

geomcols=np.arange(8)
CASSMGeom=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='Geode',header=0,usecols=geomcols)
CASSMsrc=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='CASSM_Channels',header=0,usecols=geomcols)
# # e4dhead=pd.read_excel(subterGeom,'Header')
# sourceW=pd.read_excel(subterGeom,'SourceW',header=1)
# sourceS=pd.read_excel(subterGeom,'SourceS',header=1)
# rcoordf=pd.read_excel(subterGeom,'PhonesN',header=1)
# scoordf=pd.concat([sourceW,sourceS])
# scoor=np.array(scoordf)
# rcoor=np.array(rcoordf)
# ch=rcoor[:,1].astype(int)
# ns=len(scoor[:,1])
# nr=len(rcoor[:,1])
freq=12       #frequency content (kHz)
unc=0.1    #data uncertainty (ms)
pstsrc=[6,7,8,9,10,11]
otgeom=[]
pstgeom=[]
OTrecID=[]
dfcols=['count','sourceID','receiverID','picktime','std']
OTcorr_picks=pd.DataFrame([])
PSTcorr_picks=pd.DataFrame([])

dates=['20181022085901','20181024233402','20181024235219','20181028231038',\
       '20181102232638','20181102233552','20181102235419','20181104230935',\
       '20181104234634','20181107231231','20181108003557']

OTtrigger_index=np.array([0])  
OTtrigger_ptnl_index=[] 
OTtrigger_times=np.array([0])
PSTtrigger_index=np.array([0])  
PSTtrigger_ptnl_index=[] 
PSTtrigger_times=np.array([0])

#low = 5000
#high = 8000
path=f'/home/spri902/EGS_Collab/joint_inv/'
filelist=sorted(os.listdir(f'/home/spri902/EGS_Collab/joint_inv/{dates[0]}/'))
cnt=1
#%%
for i in range(len(dates))
    for filename in sorted(os.listdir(f'/home/spri902/EGS_Collab/joint_inv/{dates[i]}/')):
        if filename.endswith(".dat"):
            sthead=_read_seg2(path+dates[i]+'/'+filename,headonly=True)
            times=sthead[0].times()
            srcCh=CASSM_chanDecode(t=times,enc=sthead[93].data)
            srcID=CASSMsrc[CASSMsrc["CytekCh"]==srcCh]
            dt=sthead[0].stats.delta 
            nyq = np.floor(0.5*(1/dt))
            st=sthead[0].stats.starttime
            et=sthead[0].stats.endtime
            sthead.trim(st,st + (dt*1999))
            OThy=sthead[12:23]
            OTtr=sthead[69:78]
            PSTtr=sthead[51:58]
            PSTtr=PSTtr + sthead[59]
            
            #OTfilt=OTtr.filter('bandpass',freqmin=low,freqmax=high,corners=4,zerophase=True)
            #PSTfilt=PSTtr.filter('bandpass',freqmin=low,freqmax=high,corners=4,zerophase=True)
            OThydata=np.stack(tr.data for tr in OThy.traces).T
            OTdataN=np.stack(tr.data for tr in OTtr.traces).T
            PSTdataN=np.stack(tr.data for tr in PSTtr.traces).T
            #sthead.detrend("linear") 
            #sthead.detrend("demean")
            #sthead.taper(type='hamming',max_percentage=0.5)
     
            slen = sthead[0].stats.npts # signal length
            t = np.arange(0, slen / (1/dt), dt) # signal time vector
            plt.figure()
            shot_gather(OTdataN,dt,trnormalize=True)
            out=plt.ginput(n=-1,timeout=0,show_clicks=True)
            for ii in range(len(out)):
                OTrecID=int(float(OTtr[ii].stats.seg2['RECEIVER_LOCATION']))
                otgeom.append((cnt,srcID.iloc[0,0],OTrecID,out[ii][1],unc))
                cnt=cnt+1

            OT_picks=pd.DataFrame(otgeom,columns=dfcols)    
            OT_picks.to_pickle(f'OT_picks_{dates[0]}')
            #OT_picks.to_csv(r'/home/spri902/EGS_Collab/joint_inv/picks/',\
                            #header=None,index=None,sep=' ',mode='a')
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            