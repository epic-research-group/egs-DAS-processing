#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 10:40:09 2021

@author: spri902
"""
import numpy as np



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





















