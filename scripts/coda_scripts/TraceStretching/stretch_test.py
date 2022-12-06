#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 07:45:21 2022

@author: spri902
"""
import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import pchip_interpolate
#%%
# this script runs an example of the moving window stretching method
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
                   
            u0tmp = u0int[startIdx : stopIdx] # locally windowed trace
            u1tmp = u1[startIdx : stopIdx] # locally windowed trace
            
            # zero lag correlation of this window  equation 6 in Hadziioannou (2012).
            ccArray[ ii, cnt ] = np.corrcoef(u0tmp.reshape(-1),u1tmp.reshape(-1))[0,1]# take cross correlation coefficient
            cnt = cnt + 1 # next window -> update counter
            
    # get the max dV within in each window
    C     = np.amax( ccArray,axis=0 )
    CCidx = ccArray.argmax(axis=0) # find the indices for those max values
    epsArray = eps[CCidx] # epsilon value at the max corrcoeff for that window
    
    return C , epsArray, tSamp
#%%
# u0 = reference trace at 6 km/s
# u1 = perturbed trace at 5.94 km/s (a -1.0% perturbation)

# load the data and plot section

data = loadmat('/home/spri902/scripts/coda_scripts/TraceStretching/traces.mat')

npts = len(data['u0']) # number of points in traces
u0 = data['u0']
u1 = data['u1']
dt = data['dt'][0][0]
tArray = np.r_[0 : npts] * dt # time vector

fig=plt.figure()
plt.plot(tArray.T,u0,label='u0')
plt.plot(tArray.T,u1,label='u1')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude [a.u.]')
plt.legend(loc='upper right')
plt.xlim([1, 3])
plt.ylim([-5e-2 ,5e-2])

# do moving window stretching

winLength = 0.5     # [s] length of moving window
tStep     = dt * 10 # [s] make a measurement every 'tStep'

dVmax     = 0.02    # set maximum for stretch parameter search (0.5 = 50%)
dV        = dVmax/6 # sample interval from -dVmax:dV:dVmax for epsilon values

ccArray, dtot, tSamp  = movingWinStretch( u0, u1, dt, winLength, tStep, dVmax, dV )

# tSamp is the array that gives the window center and changes with tStep

fig,ax = plt.subplots(2,1)
ax[0].plot( tSamp.reshape(-1), ccArray )
ax[0].set_ylabel('Corr. Coeff.') 
ax[0].set_ylim([0.95, 1])

ax[1].plot( tSamp.reshape(-1), dtot )
ax[1].set_ylabel('\epsilon')
ax[1].set_ylim([-dVmax, dVmax])
plt.xlabel('Time [s]')