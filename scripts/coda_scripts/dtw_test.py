#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 06:39:26 2022

@author: spri902
"""

import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import convolve
from coda_scripts.dtwFuncs import *

# make Ricker wavelet
def ricker(f, N, dt=0.004):
    T = dt*(N-1)
    t = np.r_[0:T+dt:dt]
    t0 = 1/f
    tau = t-t0
    # t = np.arange(-length/2, (length-dt)/2, dt)
    # y = (1.0 - 2.0*(np.pi**2)*(f**2)*(tau**2)) * np.exp(-(np.pi**2)*(f**2)*(tau**2))
    y = (1.0 - (np.pi**2)*(f**2)*(tau**2)) * np.exp(-(np.pi**2)*(f**2)*(tau**2))
    return tau, y

dt=0.004
f = 20 # A low freq of 20 Hz
# make random reflectivity sequence and convolve with Ricker
npts = 500 # length of time series
t, w = ricker(f,npts,dt)
f = np.random.rand(npts)
u0 = convolve(f,w,mode='same',method='auto')

# make time varying shifts as sine wave
amp = 0.2*dt
dp = 2*np.pi / npts
p = np.r_[0:2*np.pi:dp]
st = amp * np.sin(p) / dt

lag = 80
b=1
tvec = np.r_[0:npts]*dt
tvec2 = tvec + st
lvec = np.r_[-lag:lag+1]*dt
u1 = np.interp(tvec2,tvec,u0)
st = st/dt
stTime = st*dt

err = computeErrFunc(u1, u0, npts,lag=80)

fig1,ax = plt.subplots()

img=ax.pcolormesh(tvec,lvec,np.log10(err.T),cmap='gray')
plt.colorbar(img)
plt.xlabel('Time [s]')
plt.ylabel('Lag')
plt.title('Error Fucntion')

direction=1
dist = accumulateErrFunc(direction, err, npts,lag,b)
stbar = backtrackDistFunc(direction, dist, err, -lag, b)

stbarTime = stbar * dt
tvec2 = tvec + stbarTime

# make figure
fig2,axs = plt.subplots(2,1)
axs[0].pcolormesh(tvec,lvec,dist.T)
axs[1].plot(tvec,stTime,color='black',marker='o')
axs[1].plot(tvec,stbarTime,color='red',marker='+')
