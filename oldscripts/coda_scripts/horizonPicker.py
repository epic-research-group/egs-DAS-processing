#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 11:13:56 2022

@author: spri902
"""
import os
from datetime import datetime as datetime
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import matplotlib
%matplotlib
import matplotlib.pyplot as plt 
from matplotlib import animation
import pickle
from scipy.fft import *
from scipy.signal import hilbert, find_peaks
from operator import itemgetter
# from spectrum import dpss,pmtm


geomcols=np.arange(8)
wellcols=np.arange(12)
CASSMGeom=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='Geode',header=0,usecols=geomcols)
CASSMsrc=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='CASSM_Channels',header=0,usecols=geomcols)
CASSMsrc.columns = [c.replace(' ','_') for c in CASSMsrc.columns]
WellGeom=pd.read_excel('/home/spri902/Collab_metadata/Well_Points.xlsx',header=0,usecols=wellcols)
#%%
directory = sorted(os.walk('/home/spri902/mayCASSM'))
dates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]
xlims = mdates.date2num(dates)
stimbeg = [42, 96, 224, 352, 472, 507]
stimfin = [50,106, 232, 368, 475, 512]
stimvec=[]
[stimvec.append(np.r_[stimbeg[stim]:stimfin[stim]+1]) for stim in range(len(stimbeg))]
stimvec=np.concatenate(stimvec).ravel().tolist()
stimbeg = [96, 224, 352, 472, 507]
stimfin = [106, 232, 368, 475, 512]
stimbegLines = itemgetter(*stimbeg)(dates)
stimfinLines = itemgetter(*stimfin)(dates)
###############################################################################
# src=[0] # OB source
# # OT hydrophone and accels listed by channel number ** not index **
# recVec = np.array((13, 14, 15, 16, 17, 18, 19, 20, 21,
#           22, 23, 70, 71, 72, 73, 74, 75, 76, 77, 78))
# # channels for OT hydrophones and accels
# hydChan=slice(12,23)
# accChan=slice(69,78)
###############################################################################
# PST src and channels for PDB hydrophones and accels
# src=[17] # PST source
# hydChan=slice(0,12) # PDB hydrophones
# accChan=slice(30,36)# PDB accels
# ac2Chan=slice(39,42)# PDB accels
src=[0] # source number
accChan = slice(54,57)
recVec = np.array((55,56,57))
# hydChan=slice(0,12) # PDB hydrophones
# accChan=slice(24,30)# PDB accels
# ac2Chan=slice(36,39)# PDB accels
# PDB hydrophone and accels listed by channel number ** not index **
# recVec = np.array((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
#                     31, 32, 33, 34, 35, 36, 40, 41, 42))
# recVec = np.array((25,26,27,28,29,30,37,38,39))

###############################################################################
# pull src and rec locations
srcloc = CASSMsrc.loc[src[0],'Easting_m':'Elev_m'].values
srcloc = srcloc.reshape(-1,3).astype(float)
recloc = CASSMGeom.loc[recVec-1,'Easting(m)':'Elev(m)'].values
recInd = np.arange(len(recVec))
# only need to grab part of the signal (don't need all 9600 samples)
siglen = 2000
samples = np.arange(siglen)

#%% Load data
filepath = '/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/src0_PSTacc/'
datfile = 'OB0data.npy'
seisdata = np.load(os.path.join(filepath, datfile))
os.chdir(filepath)
# for i in recInd:
#     recID = CASSMGeom[CASSMGeom["Channel"] == recVec[i]].Sensor.item()
#     fig = pickle.load(open(f'CASSM_Source_{src[0]}_Receiver_{recID}','rb'))
xcorrA = []
nx=seisdata.shape[1]
lags = np.arange(-nx + 1, nx)
for j in range(seisdata.shape[2]):
    maxlag=[]
    ccorr=[]
    recID = CASSMGeom[CASSMGeom["Channel"] == recVec[j]].Sensor.item()
    for i in range(1,seisdata.shape[0]):
        xcorr = np.correlate(seisdata[0,:,j], seisdata[i,:,j],mode='full')
        xcorr /= xcorr[nx - 1]
        maxlag.append(lags[np.argmax(xcorr)])
        ccorr.append(np.corrcoef(seisdata[0,:,j],seisdata[i,:,j])[0,1])
        # if i%20 == 0:
        #     fig,ax = plt.subplots()
        #     ax.plot(lags,xcorr,'r')
        #     ax.set_xlabel('lag')
        #     ax.set_ylabel('corr coeff')
        #     ax.grid(True)
    xcorrA.append(ccorr)        
    fig2,ax = plt.subplots()
    ax.set_xlim(dates[0],dates[-1])
    plt.plot(dates[0:len(ccorr)],ccorr,'.k')
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
    ax.xaxis_date()
    plt.gcf().autofmt_xdate()
    ax.set_ylim(-1,1)
    ax.set_xlabel('Dates')
    ax.set_ylabel('Normalized CorrCoeff')
    plt.title(f'Zero-lag CorrCoef CASSM Source: {src[0]} Receiver_{recID}')
    pickle.dump(ax, open(f'Zero-lag_CorrCoef_CASSM-Source_{src[0]}_Receiver_{recID}','wb'))
    plt.close()
#%%
fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(12,6))
plt.suptitle('CASSM Accelerometer Cross-Correlations')
a1 = ax1.plot(dates[:-1],ccorr[4,:],'.k',label='Hydro #4 Xcorr')
ax1.set_ylabel('Zero-lag Xcorr')
l4 = ax1.axvspan(stimbegLines[0],stimfinLines[0],alpha=0.7,color='grey',label='Injection')
l5 = ax1.axvspan(stimbegLines[1],stimfinLines[1],alpha=0.7,color='grey',label='Injection')
l6 = ax1.axvspan(stimbegLines[2],stimfinLines[2],alpha=0.7,color='grey',label='Injection')
l7 = ax1.axvspan(stimbegLines[3],stimfinLines[3],alpha=0.7,color='grey',label='Injection')
l8 = ax1.axvspan(stimbegLines[4],stimfinLines[4],alpha=0.7,color='grey',label='Injection')
lns2 = a1+[l5]
lbls2 = [ln.get_label() for ln in lns2]
ax1.legend(lns2,lbls2,loc='lower left')

a2 = ax2.plot(dates[:-1],ccorr[15,:],'.k',label='Accel Xcorr')
ax2.set_ylabel('Zero-lag Xcorr')
l4 = ax2.axvspan(stimbegLines[0],stimfinLines[0],alpha=0.7,color='grey',label='Injection')
l5 = ax2.axvspan(stimbegLines[1],stimfinLines[1],alpha=0.7,color='grey',label='Injection')
l6 = ax2.axvspan(stimbegLines[2],stimfinLines[2],alpha=0.7,color='grey',label='Injection')
l7 = ax2.axvspan(stimbegLines[3],stimfinLines[3],alpha=0.7,color='grey',label='Injection')
l8 = ax2.axvspan(stimbegLines[4],stimfinLines[4],alpha=0.7,color='grey',label='Injection')
lns2 = a2+[l5]
lbls2 = [ln.get_label() for ln in lns2]
ax2.legend(lns2,lbls2,loc='lower left')

a3 = ax3.plot(dates[:-1],ccorr[11,:],'.k',label='Hydro #12 Xcorr')
ax3.set_ylabel('Zero-lag Xcorr')
ax3.set_xlabel('Date')
l4 = ax3.axvspan(stimbegLines[0],stimfinLines[0],alpha=0.7,color='grey',label='Injection')
l5 = ax3.axvspan(stimbegLines[1],stimfinLines[1],alpha=0.7,color='grey',label='Injection')
l6 = ax3.axvspan(stimbegLines[2],stimfinLines[2],alpha=0.7,color='grey',label='Injection')
l7 = ax3.axvspan(stimbegLines[3],stimfinLines[3],alpha=0.7,color='grey',label='Injection')
l8 = ax3.axvspan(stimbegLines[4],stimfinLines[4],alpha=0.7,color='grey',label='Injection')
lns2 = a3+[l5]
lbls2 = [ln.get_label() for ln in lns2]
ax3.legend(lns2,lbls2,loc='lower left')
#%% Save zero lag ccorr data
with open('OB0_ccorr.npy','wb') as f:
    np.savez(f,xcorrA)
    
ccorrdata = np.load('PST9_ccorr.npy')
ccorr = ccorrdata['arr_0']

#%%     
for i in range(seisdata.shape[2]):
    recID = CASSMGeom[CASSMGeom["Channel"] == recVec[i]].Sensor.item()
    fig = pickle.load(open(f'Zero-lag_CorrCoef_CASSM-Source_{src[0]}_Receiver_{recID}','rb'))
#%% Make instantaneous amplitude, phase and frequency calcs and plots
# datslice = slice(195,1351)
datrange = np.r_[0:2000]
sigE     = []
instF = []
phsyn    = []

fs       = 48000
dt=1/fs
# time half bandwidth parameter for pmtm
# NW = 3.5
# # uses the first k number of slepian sequences
# k  = 4
# [tapers, eigen] = dpss(siglen, NW, k)
# Sk, weights, eigvals=pmtm(recArray[:,0], e=eigen, v=tapers, NFFT=2048, show=False)
# Sx = abs(Sk)**2
# Sx = np.mean(Sx * np.transpose(weights), axis=0) * dt

# yf = fft(recArray)
# xf = fftfreq(siglen,dt)

for i in range(len(recVec)):

    recArray = seisdata[:,:,i].T
    recArray = recArray - recArray.mean(axis=0)

    a = hilbert(recArray) # make the analytic signal for the traces
    # the Hilbert transformed signal can be obtained from np.imag(hilbert(x))
    # the original signal can be obtained from np.real(hilbert(x))

    
    # e = np.abs(a) # get the magnitude of the analytsig for amplitude envelope

    p = np.unwrap(np.angle(a,deg=False)) # get instantaneous phase
    # pp = np.angle(a,deg=False) # unwrapped instantaneous phase

    f = (np.diff(p,prepend=0) / (2.0*np.pi) * fs) # get instantaneous frequency
    T = 1/f
    T = np.diff(T,prepend=0)
    # refT = recArray[:,j-1]
    # warT = recArray[:,j]
    # phase_sync = 1 - np.sin(np.abs(p-p[:,0,None]) / 2)
    # phsyn.append(phase_sync)
    # sigE.append(e)
    # instF.append(f)
    
    # make a ref trace to warped trace plot
    # fig,(ax0,ax1,ax2,ax3,ax4,ax5) = plt.subplots(6,1,figsize=(15,10),sharex=True)
    # ax0.plot(refT,color='r',label='ref')
    # ax0.plot(warT,color='b',label='warped')
    # ax0.set(title='Timeseries Data')
    # ax0.legend()
    
    # ax1.plot(refT,color='r',label='ref')
    # ax1.plot(e1,color='g',label='ref envelope')
    # ax1.legend()
    # ax2.plot(warT,color='b',label='warped')
    # ax2.plot(e2,color='m',label='warped envelope')
    # ax2.legend()
    
    # ax3.plot(p1,color='r',label='ref')
    # ax3.plot(p2,color='b',label='warped')
    # ax3.set(title='Instantaneous Phase')
    # ax3.legend()
    
    # ax4.plot(f1,color='r',label='ref')
    # ax4.plot(f2,color='b',label='warped')
    # ax4.set(title='Instantaneous frequency')
    # ax4.legend()
    
    # ax5.plot(phase_sync)
    # ax5.set(ylim=[0,1.1],title='Ref-Warped Phase Synchrony')
    # plt.tight_layout()
    # ax5.legend()
        
    # make a cumulative plot comparing all warped traces to ref trace
    fig,ax = plt.subplots()
    recID = CASSMGeom[CASSMGeom["Channel"] == recVec[i]].Sensor.item()
    vm=np.percentile(T,99)
    img=ax.pcolormesh(xlims,datrange,T,cmap='RdBu',vmin=-vm,vmax=vm)
    # img=ax.pcolormesh(np.arange(1,seisdata.shape[0]),datrange,np.array(sigE).T,cmap='RdBu',vmin=-vm,vmax=vm)
    # img=ax.pcolormesh(xlims,datrange,np.array(instFreq).T,cmap='RdBu',vmin=-vm,vmax=vm)
    plt.gca().invert_yaxis()
    plt.colorbar(img)
    plt.title(f'PhaseTracking_CASSM_Source_{src[0]}_Receiver_{recID}')
    ax.set_xlim(dates[0],dates[-1])
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
    ax.xaxis_date()
    plt.gcf().autofmt_xdate()
    
    phase_sync = []
    e =[]
    f=[]
    # pickle.dump(ax, open(f'PhaseTracking_CASSM-Source_{src[0]}_Receiver_{recID}','wb'))
    # plt.close()
    
    
for i in range(seisdata.shape[2]):
    recID = CASSMGeom[CASSMGeom["Channel"] == recVec[i]].Sensor.item()
    fig = pickle.load(open(f'PhaseTracking_CASSM-Source_{src[0]}_Receiver_{recID}','rb'))
#%% load stretching and xcorr data
stretchdata = np.load('PST11_stretch_16.npz')
dtArray = stretchdata['arr_0']
tSamp = stretchdata['arr_1']
cArray = stretchdata['arr_2']

ccorrdata = np.load('PST11_ccorr.npy')
ccorr = ccorrdata['arr_0']
t= np.r_[0:2000]*dt
#%% Plot single receiver waveforms through time
# animation line plot example

fig, ax = plt.subplots(2, 1, figsize = (16, 6))
ax[1] = ax[0].twinx()
def animate(i):
    
    ax[0].cla() # clear the previous image
    # ax[1].cla() # clear the previous image
    
    if i in stimvec:
        ax[0].plot(t,seisdata[0,:,3],label='Reference Trace')
        ax[0].plot(t,seisdata[i,:,3],color='c',label='Perturbed Trace',alpha=0.8) # plot the line
        ax[0].set_xlim([t[0], t[-1]]) # fix the x axis
        plt.title(f'CASSM Data Acquisition Time {dates[i]}')
        plt.legend()
    else:
        ax[0].plot(t,seisdata[0,:,3],label='Reference Trace')
        ax[0].plot(t,seisdata[i,:,3],label='Perturbed Trace',alpha=0.8) # plot the line
        ax[0].set_xlim([t[0], t[-1]]) # fix the x axis
        plt.title(f'CASSM Data Acquisition Time {dates[i]}')
        plt.legend()
    
    # ax[1].plot(ccorr[17,i],'.k')
    # ax[1].plot(tSamp,dtArray[i,:])
    ax[1].plot(tSamp,dtArray)
    # ax[1].set_xlim([0, tSamp[-1]]) # fix the x axis
    
    
anim = animation.FuncAnimation(fig, animate, frames = seisdata.shape[0], interval = 1, repeat=False,blit = False)
plt.show()


