#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 14:45:15 2022

%
% USAGE: [C,epsArray,tSamp] = movingWinStretch(u0,u1,dt,twin,dVmax,dV)
% INPUT:
%   u0    = reference trace
%   u1    = perturbed trace
%   dt    = sample interval [s]
%   twin  = length of analysis window [s]
%   tstep = time step of windows [s]
%   dVmax = maximum epsilon value to test (i.e. 0.5 = 50%)
%   dV    = sample interval for epsilon vector (i.e. epsilon = -dVmax:dV:dVmax)
% OUTPUT:
%   C        = maximum correlation coefficient in each window
%   epsArray = epsilon value for each window corresponding to maximum correlation coefficient in C
%   tSamp    = time vector containing the center time of each analyzed window
%
% This function computes the moving window stretching method to estimate
% velocity changes.
%
% 
% Written by Dylan Mikesell (mikesell@mit.edu)
% Last modified 13 April 2015
%
% Guts of algorithm taken from an example by Berenice Froment
% Ported to Python by Parker Sprinkle March 2018
"""
import matplotlib
%matplotlib
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
from scipy.interpolate import pchip_interpolate
import os
import matplotlib.dates as mdates
from datetime import datetime as datetime
import scipy.io
#%% Read in stim flow data
# sf0523 = pd.read_csv('~/EGS_Collab/stimflow/SF_E1_Stim2_1.5mStepTest.csv',
#                      header=[0,1],dtype={'Quizix Flow':np.float64})
# sf0524 = pd.read_csv('~/EGS_Collab/stimflow/SF_E1_Stim2_DriveToProduction.csv',
#                      header=[0,1],dtype={'Triplex Net Flow':np.float64})
# sftime = pd.concat([sf0523.Time, sf0524.Time])
# sfdat = pd.concat([sf0523.iloc[:,2],sf0524.iloc[:,4]])


#%%
def movingWinStretch(u0,u1,dt,twin,tstep,dVmax,dV):

    npts  = len( u1 )         # number of data points in the traces
    nstep = np.round( tstep / dt)  # number of samples to step through windows
    nwin  = np.round( twin / dt )  # number of samples in the window
    nwin2 = np.floor( nwin / 2 )   # half the sample window (in samples)
    eps   = np.arange(-dVmax, dVmax+dV, dV) # dV array to test
    neps  = len( eps )        # number of dV values to test

    nMeas = len( np.arange(nwin2 + 1, (npts - nwin2)+1, nstep) ) # number of measurements

    ccArray = np.zeros((neps, nMeas)) # allocate correlation coefficient vector
    tSamp   = np.zeros((1, nMeas)) # allocate window center time vector
    
    for ii in range(neps): # test each epsilon value
        
        print( f'Testing epsilon = {eps[ii]}' )
        
        # stretch the reference trace for the given eps_vec
        tPrime =  np.r_[0 : npts] * ( 1 + eps[ii] ) # stretch the time axis
        u0int  = pchip_interpolate(np.r_[0:npts],u0,tPrime,axis=0) # interpolate the stretched trace
        cnt = 0 # counter for windows
        
        # jj is the sample number at the window center
        for jj in np.arange(nwin2 + 1, npts - nwin2, nstep): # loop through windows 

            
            tSamp[0][cnt] = (jj-1) * dt # [s] center of window  
            
            startIdx = int(jj - nwin2 - 1) # window starts half width before center of window
            stopIdx  = int(jj + nwin2 - 1) # window ends half width after center of window
                   
            u0tmp = u0int[ startIdx : stopIdx ] # locally windowed trace
            u1tmp = u1[    startIdx : stopIdx ] # locally windowed trace
            
            # zero lag correlation of this window  equation 6 in Hadziioannou (2012).
            ccArray[ ii, cnt ] = np.corrcoef(u0tmp.reshape(-1),u1tmp.reshape(-1))[0,1]# take cross correlation coefficient
            cnt = cnt + 1 # next window -> update counter
            
    # get the max dV within in each window
    C     = np.amax( ccArray,axis=0 )
    CCidx = ccArray.argmax(axis=0) # find the indices for those max values
    epsArray = eps[CCidx] # epsilon value at the max corrcoeff for that window
    
    return C , epsArray, tSamp
#%%
# traces = scipy.io.loadmat('/home/spri902/EGS_Collab/traces.mat')
# dt=traces['dt']
# u0=traces['u0']
# u1=traces['u1']
# winlen=0.5
# tStep=dt*10
# dVmax = 0.02
# dV = dVmax/6
# ccArray, dtot, tSamp  = movingWinStretch( u0, u1, dt, winlen, tStep, dVmax, dV )

filepath = '/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/src0_PSTacc/'
datfile = 'OB0data.npy'
seisdata = np.load(os.path.join(filepath,datfile))
directory = sorted(os.walk('/home/spri902/mayCASSM'))
dates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]
xlims = mdates.date2num(dates)

dt    = 1/48000
twin  = dt*80
tstep = dt*3
dVmax = 0.5
dV    = dVmax/12
t = np.r_[0:seisdata.shape[1]]*dt
cArray = []
dtArray = []
tArray = []
for i in range(1,seisdata.shape[0]):
    u0 = seisdata[0,:,2] # use 3 or 16 depending on which accel (PDT or PDB)
    u1 = seisdata[i,:,2]
    
    C, dtot, tSamp = movingWinStretch(u0, u1, dt, twin, tstep, dVmax, dV)
    cArray.append(C)
    dtArray.append(dtot)
    tArray.append(tSamp)
    
    # fig = plt.figure(figsize=(10,5))
    # ax0 = fig.add_subplot(3,1,1)
    # ax1 = fig.add_subplot(3,1,2, sharex = ax0)
    # ax2 = fig.add_subplot(3,1,3, sharex = ax1)
    # ax0.plot(t,u0,t,u1)
    # ax1.plot(tSamp.reshape(-1),C)
    # ax1.set(title='Corr. Coeff')
    # ax2.plot(tSamp.reshape(-1),dtot)
    # ax2.set(xlim=[0,dt*npts],ylim=[-dVmax,dVmax],title='Epsilon : dt/t')
    # plt.tight_layout()
    # plt.xlabel('Time (s) ')

tSamp = tSamp.reshape(-1)
fig, ax = plt.subplots()
# Twin the x-axis twice to make independent y-axes.
axes = [ax, ax.twinx()]

# Make some space on the right side for the extra y-axis.
# fig.subplots_adjust(right=0.75)
# Move the last y-axis spine over to the right by 20% of the width of the axes
# axes[-1].spines['right'].set_position(('axes', 1.2))

# To make the border of the right-most axis visible, we need to turn the frame
# on. This hides the other plots, however, so we need to turn its fill off.
# axes[-1].set_frame_on(True)
# axes[-1].patch.set_visible(False)
img=axes[0].pcolormesh(xlims[1:],tSamp,np.array(dtArray).T*100,cmap='RdBu')
plt.gca().invert_yaxis
axes[0].xaxis_date()
plt.gcf().autofmt_xdate()
plt.title(f'Trace Streching for Source {datfile}')
axes[0].set_xlabel('Experiment Date')
axes[0].set_ylabel('Trace Time')
axes[0].set_ylim([0.039,0.0025])

# l1=axes[1].plot(xlims[1:],ccorr[16,:],color='black')
# l2=axes[2].plot(dasdnums,PDBdat[29,:],color='black')
axes[1].set_ylabel('Correlation Coefficients')

cbar=plt.colorbar(img, orientation='horizontal',pad=0.2)
cbar.set_label('dV -> Velocity Change (%)')


#%% Save and load data
os.chdir(filepath)
with open('PST9_stretch_16.npz','wb') as f:
    np.savez(f,dtArray,tSamp,cArray)
    
stretchdata = np.load('PST9_stretch_16.npz')
dtArray = stretchdata['arr_0']
tSamp = stretchdata['arr_1']
cArray = stretchdata['arr_2']

ccorrdata = np.load('PST9_ccorr.npy')
ccorr = ccorrdata['arr_0']
