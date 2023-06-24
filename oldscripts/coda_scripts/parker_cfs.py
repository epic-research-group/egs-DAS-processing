#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:33:40 2020

@author: spri902
"""



# from obspy import read
import numpy as np
from scipy.stats import kurtosis
from scipy.spatial.distance import cdist
import scipy.signal as signal


def aic_maeda(dat):
    num_samples = dat.shape[0]
    aic = np.zeros_like(dat)
    for idx in range(num_samples):
        aic[idx] = idx*np.log(dat[:idx+1].var()) + \
            (num_samples-idx)*np.log(dat[idx:].var())
    aic[-1] = aic[-2]
    return aic



def kurtosis_fisher(dat):
    return (((dat-dat.mean())/dat.std())**4).mean() - 3

def compute_kurtosis(dat, win):   
    kurt = np.zeros_like(dat)
    gkurt= np.zeros_like(dat)
    for idx in range(win, kurt.size):
        kurt[idx] = kurtosis_fisher(dat[idx-win:idx]) 
        
    gkurt=np.gradient(kurt)                 
    return gkurt

def aic_parker(dat):
    npts = dat.shape[0]
    AIC = np.zeros(npts)
    
    # reverse indexing to remove the nan, np.std(data[:0]) is nan, starting index need to 
    # be npts-2, if data array only has 1 sample, the std is 0, the log10(0) is inf
    for k in range(npts-2,0,-1):

      a = k*np.log10(np.std(dat[:k])**2)+(npts-k-1)*np.log10(np.std(dat[k:])**2)
      
      #print a,np.log10(np.std(data[k:]))
      if a == -float('inf'):
          a = AIC[k+1]
      AIC[k] = a
      AIC[0] = AIC[1]
      AIC[-1] = AIC[-2]
      AIC = np.array(AIC)
    
      AIC_deriv = []
    for i in range(npts-1):
      b = np.abs(AIC[i+1]-AIC[i])
      AIC_deriv.append(b)

    AIC_deriv.insert(0,AIC_deriv[0])
    AIC_deriv = np.array(AIC_deriv)

    return AIC_deriv


def signaltonoise(a, axis=0, ddof=0):
    a = np.asanyarray(a)
    m = a.mean(axis)
    sd = a.std(axis=axis, ddof=ddof)
    return np.where(sd == 0, 0, m/sd)

def threshold(trace,t_ma,nsigma):
    """ 
    Control the threshold level with nsigma. Pulled from Austin Holland's PhasePapy code'
    """
    dt = trace.stats.delta
    npts_Tma = int(round(t_ma / dt,0))
    LEN = trace.stats.npts
    

    threshold = np.zeros(LEN)
    threshold[npts_Tma:LEN] = rms(rolling_window(np.array(aic_parker(trace.data[0:LEN-1])).reshape(-1),npts_Tma), -1) * nsigma
    #threshold[npts_Tma:LEN] = rms(rolling_window(np.array(compute_kurtosis(trace.data[0:LEN-1],5)).reshape(-1),npts_Tma), -1) * nsigma
    threshold[0:npts_Tma] = 1000

    return threshold

def rms(x, axis=None):
    """ Function to calculate the root mean square value of an array.
    """
    return np.sqrt(np.mean(x**2, axis=axis))   
    
def rolling_window(a, window):
    """ Efficient rolling statistics with NumPy: This is applied to 
        threshold() to calcuate threshold to locate signal arrival
        Reference from:
        http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides) 


def reject_outliers(data, m = 4.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

def set_bit(v, index, x):
  """Set the index:th bit of v to 1 if x is truthy, else to 0, and return the new value."""
  mask = 1 << index   # Compute mask, an integer with just bit 'index' set.
  v &= ~mask          # Clear the bit indicated by the mask (if x is False)
  if x:
    v |= mask         # If x was True, set the bit indicated by the mask.
  return v            # Return the result, we're done.

def CASSM_chanDecode(t,enc,fig=None):
    """
% Notes
% -------------------
%
% Decodes binary CASSM encoding for source identification. 
%
% Input Arguments
% -------------------
% 1. t         :  array of time values for the encoding data (in ms (!))
% 2. enc       :  the CASSM encoding channel data
% 3. threshold :  clips data below threshold (use .5)
% 4. fig       :  plots results on figure(fig) id fig>0
%  
% Output Variables
% -------------------
% 1. channel   :  the decoded Cytek channel number
%
% History
% -------------------
% 1. First version by Todd Wood
% 2. Edits by JBAF : April 26th, 2018
% 
% Modifications
% -------------------
% 2. Clean-up and documentation during debugging process. Took out a
% "round" statement which was always pushing bit_width to zero (right
% move)? Might need to run this by Todd.
%
%
    """

    BYTE_SIZE=5
    bit_threshold=0.5
    enc[np.argwhere(enc < bit_threshold*np.max(enc))] = 0
    enc[np.argwhere(enc>0)] = np.max(enc)

    # calculating 1st derivative to get bit boundaries
    fDiff = np.abs(np.diff(enc))

    #finding the times for locations with a high derivative
    b = t[np.argwhere(fDiff > bit_threshold*np.max(fDiff))]

    # calculating the width and start location from the first peak
    # bit_width = np.round(b[1] - b[0])
    bit_width = (b[1] - b[0]) #bit_width in time
    bit_start = b[0]

    # this must assume a guard bit - spaced afterwards ....
    data_start = bit_start + (1.5*bit_width)
    bits = []
    channel = 0

    # looping over "bits" in signal
    for k in range(BYTE_SIZE):
        # finding bit location
        bits.append(data_start + (k+1)*bit_width)
        if bits[k] > t[-1]:
            bits[k] = t[-1]
        if enc[np.argwhere(np.min(t[t >= bits[k]])==t)] > 0:
            channel = set_bit(channel,k,1)
    return channel        

    # simple plot to verify
    if fig:
        plt.figure(fig)
        plt.plot(t,enc)
        axes=plt.gca()
        axes.set_xlim([np.ceil(data_start-2*bit_width), data_start+8*bit_width])
        plt.plot(t,np.append(np.abs(np.diff(enc)),0),'g-')
        plt.plot(bits,np.ones((1,len(bits))).reshape((5,1)) * bit_threshold*np.max(enc),'r+')
        plt.xlabel('Time')
        plt.ylabel('Signal')
        
def kurt_parker(dat,dt):
    m = len(dat)
    fs=1/dt
    t_win=t_ma=(m)*(dt)/1500
    t = np.arange(0, dt * m, dt)
    
    Nsta = int(t_win * fs)

    # compute the short time average (STA)
    kt = np.zeros(m, dtype='float64')
    pad_kt = np.zeros(Nsta)
    # Tricky: Construct a big window of length len(a)-nsta. Now move this
    # window nsta points, i.e. the window "sees" every point in a at least
    # once.
    # Changed xrange to range as it is compatible in both python 2 & 3
    for i in range(m):  # window size to smooth over
        kt[i] = abs(kurtosis(dat[i-Nsta:i]))

    kt[0:Nsta] = 0
    
    return kt
        
def sl_velocity(src_locs,rec_locs,arv_times):    
    dist = cdist(src_locs,rec_locs,'euclidean')
    
    vel = dist/arv_times
    return dist , vel    
def find_outlierIQR(df):
    Q1 = df.quantile(0.25).picktime
    Q3 = df.quantile(0.75).picktime
    IQR = Q3-Q1
    ol = df.loc[(df.picktime < (Q1-1.5*IQR)) | (df.picktime > (Q3+1.5*IQR))]
    return ol     
def find_outlierHampel(df):
    med = df.median().picktime
    nl = abs(df.picktime - med)
    cond = nl.median()*4.5
    gl = df.loc[(nl>cond)]
    return gl


        
        
        
        
        
        
        
        
        
        
        