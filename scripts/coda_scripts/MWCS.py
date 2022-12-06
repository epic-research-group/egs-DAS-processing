#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon September 12, 2022

%
% USAGE: [C,epsArray,tSamp] = MWCS(u0,u1,dt,twin,dVmax,dV)
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
% This function computes the moving window cross spectral (MWCS) method to estimate
% velocity changes.
%
% Written by D. Parker Sprinkle 
"""


import matplotlib
%matplotlib
import matplotlib.pyplot as plt 
import numpy as np
from obspy.signal.invsim import cosine_taper
from obspy.signal.regression import linear_regression
from scipy.signal import  detrend, fft, fftfreq, boxcar, hanning
import os
import matplotlib.dates as mdates
from datetime import datetime as datetime

#%% Function defs
def nextpow2(x):
    """
    Returns the next power of 2 of x
    """
    return np.ceil(np.log2(np.abs(x)))

def smooth(x, window='boxcar', half_win=3):
    """ some window smoothing """
    window_len = 2 * half_win + 1
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    if window == "boxcar":
        w = boxcar(window_len).astype('complex')
    else:
        w = hanning(window_len).astype('complex')
    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[half_win:len(y) - half_win]

def getCoherence(dcross, ref, curr):
    n = len(dcross)
    coh = np.zeros(n).astype('complex')
    valids = np.argwhere(np.logical_and(np.abs(ref) > 0, np.abs(curr > 0)))
    coh[valids] = dcross[valids] / (ref[valids] * curr[valids])
    coh[coh > (1.0 + 0j)] = 1.0 + 0j
    return coh


def movingWinCrossSpec(u0,u1,freqmin,freqmax,fs,tmin,winlen,tstep,smoothing_win=5):
    '''
:type u0: :class:`numpy.ndarray`
:param u0: The "Current" timeseries
:type u1: :class:`numpy.ndarray`
:param u1: The "Reference" timeseries
:type freqmin: float
:param freqmin: The lower frequency bound to compute the dephasing (in Hz)
:type freqmax: float
:param freqmax: The higher frequency bound to compute the dephasing (in Hz)
:type fs: float
:param fs: The sampling rate of the input timeseries (in Hz)
:type tmin: float
:param tmin: The leftmost time lag (used to compute the "time lags array")
:type winen: float
:param winlen: The moving window length (in seconds)
:type tstep: float
:param tstep: The step to jump for the moving window (in seconds)
:type smoothing_win: int
:param smoothing_win: If different from 0, defines the half length of
    the smoothing hanning window.
    '''

    delta_t    = []
    delta_err  = []
    delta_mcoh = []
    time_axis  = []
    
    window_len_samples = np.int(winlen * fs)
    padz = int(2** (nextpow2(winlen) + 2))

    cnt =0
    coswin = cosine_taper(winlen,0.85)
    minInd = 0
    maxInd = winlen

    while maxInd <= len(u1):
        u1tmp  = u1[minInd:(minInd + winlen)]
        u1tmp  = detrend(u1tmp,type='linear')
        u1tmp *= coswin

        u0tmp  = u0[minInd:(minInd + winlen)]
        u0tmp  = detrend(u0tmp,type='linear')
        u0tmp *= coswin

        minInd += int(tstep * fs)
        maxInd += int(tstep * fs)

        u0fft = fft(u0tmp,n=padz)[:padz // 2]
        u1fft = fft(u1tmp,n=padz)[:padz // 2]

        u0f2 = np.real(u0fft) ** 2 + np.imag(u0fft) ** 2
        u1f2 = np.real(u1fft) ** 2 + np.imag(u1fft) ** 2

        # compute the cross spectrum
        Xov = np.multiply(u0,np.conjugate(u1))
        if smoothing_win != 0:
            du1 = np.sqrt(smooth(u1f2,window='hanning',
                                half_win=smoothing_win)) 
            du0 = np.sqrt(smooth(u0f2,window='hanning',
                                half_win=smoothing_win))
            X = smooth(X,window='hanning',
                            half_win=smoothing_win)

        else:
            du1 = np.sqrt(u1f2)
            du0 = np.sqrt(u0f2)

        dcross = np.abs(X)

        # find values in cross spectrum in the right freq bands
        freq_vec    = fftfreq(len(X) * 2, 1. / fs)[:padz // 2]
        index_range = np.argwhere(np.logical_and(freq_vec >= freqmin,
                                                freq_vec <= freqmax)) 

        # Now find the coherence
        coh     = getCoherence(dcross,du0,du1)
        corange = coh[index_range]
        meancoh = np.mean(corange)

        # get weights
        w                  = 1.0 / (1.0 / (corange ** 2) - 1.0)
        w[corange >= 0.99] = 1.0 / (1.0 / 0.98101 - 1.0)
        w                  = np.sqrt(w * np.sqrt(dcross[index_range]))
        w                  = np.real(w)

        # Freq array
        v = np.real(freq_vec[index_range]) * 2 * np.pi

        # phase
        phi    = np.angle(X)
        phi[0] = 0.
        phi    = np.unwrap(phi)
        phi    = phi[index_range]

        # Calculate the slope with a weighted least square linear regression
        # forced through the origin
        # weights for the WLS must be the variance !
        m, em = linear_regression(v.flatten(), phi.flatten(), w.flatten())

        delta_t.append(m)

        # print phi.shape, v.shape, w.shape
        e = np.sum((phi - m * v) ** 2) / (np.size(v) - 1)
        s2x2 = np.sum(v ** 2 * w ** 2)
        sx2 = np.sum(w * v ** 2)
        e = np.sqrt(e * s2x2 / sx2 ** 2)

        delta_err.append(e)
        delta_mcoh.append(np.real(meancoh))
        time_axis.append(tmin + winlen / 2. + count * tstep)
        count += 1

        del u0fft,u1fft
        del X 
        del freq_vec
        del index_range
        del w, v, e, s2x2, sx2, m, em
    if maxInd > len(u1) + tstep * fs:
        print("The last window was too small but was computed")

    return np.array([time_axis,delta_t,delta_err,delta_mcoh]).T



#%% Setup variables for MWCS analysis
filepath = '/home/spri902/EGS_Collab/4850/results/maystim/processed_CASSM/single_src_rec_gathers/src9_PDBhyds_accs/'
datfile = 'PST9data.npy'
seisdata = np.load(os.path.join(filepath,datfile))
directory = sorted(os.walk('/home/spri902/mayCASSM'))
dates = [datetime.strptime(d,'%Y%m%d%H%M%S') for d in sorted(directory[0][1])]
xlims = mdates.date2num(dates)

dt    = 1/48000 # sampling time interval
twin  = dt*80 # window length in samples (num of samples in window)
tstep = dt*3 # time steps to slide the moving window
dVmax = 0.5 
dV    = dVmax/12
t = np.r_[0:seisdata.shape[1]]*dt
cArray = []
dtArray = []
tArray = []
