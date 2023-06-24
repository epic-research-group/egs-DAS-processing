#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:46:41 2020

@author: spri902
"""
# Use the preceding window to the peak of stacked trace criterion based on signal observations. 
# Since all EMP have higher frequencies and arrive at the same time on all receivers, 
# try applying a STA/LTA or kurtosis on a high-pass filtered waveform envelope stack to determine the onset of EMP.
import numpy as np
import obspy
from obspy.io.seg2.seg2 import _read_seg2
from scipy.signal import *
from phasepickers import aic_parker

# def empPick(stream):
#     times=stream[0].times()

#     stream=stream.detrend('demean')
#     stream=stream.detrend('linear')
#     empcopy=stream.copy()

#     empstack=empcopy.stack(group_by='all',stack_type='linear',npts_tol=1,time_tol=0)

#     # win=50
#     # highkurt, gradkurt=compute_kurtosis(empstack[0].data[4800:4900], win)
#     # gradpeaks, _ =find_peaks(gradkurt[4790:4890],height=4)[:]
#     # kurtpeaks, _=find_peaks(highkurt,height=30)[:]
#     empaic,empaicd=aic_parker(empstack[0].data[4798:4898])
#     empPickind=4798+np.argmax(empaicd)-1

#     empLB=4804
#     empUB=empLB + 80
        
#     # peakrefine=empLB+np.argmax(np.abs(empstack[0].data[empLB:empUB]))
#     peaks=find_peaks(empstack[0].data[empLB:empUB])
#     proms=peak_prominences(empstack[0].data[empLB:empUB],peaks[0])
#     peakrefine=empLB + peaks[0][np.argmax(proms[0])]
#     return peakrefine, times[peakrefine],empstack




def empPick2020(stream):
    times=stream[0].times()

    stream=stream.detrend('demean')
    stream=stream.detrend('linear')
    empcopy=stream.copy()

    empstack=empcopy.stack(group_by='all',stack_type='linear',npts_tol=1,time_tol=0)

    # win=50
    # highkurt, gradkurt=compute_kurtosis(empstack[0].data[4800:4900], win)
    # gradpeaks, _ =find_peaks(gradkurt[4790:4890],height=4)[:]
    # kurtpeaks, _=find_peaks(highkurt,height=30)[:]
    empaic,empaicd=aic_parker(empstack[0].data[480:580])
    empPickind=480+np.argmax(empaicd)-1

    empLB=480
    empUB=empLB + 85
        
    # peakrefine=empLB+np.argmax(np.abs(empstack[0].data[empLB:empUB]))
    peaks=find_peaks(empstack[0].data[empLB:empUB])
    proms=peak_prominences(empstack[0].data[empLB:empUB],peaks[0])
    peakrefine=empLB + peaks[0][np.argmax(proms[0])]
    return peakrefine, times[peakrefine],empstack, empPickind,times[empPickind]





# data = np.stack(tr.data for tr in sthead.traces)
# data=data.T
# filter data
# filtcopy=sthead.copy()
# bpfilttraces=[obspy.signal.filter.bandpass(filtcopy[ii].data,freqmin=100,freqmax=2500,df=48000,corners=5, zerophase=True) for ii in range(0,len(sthead))]

# hifiltertraces=[obspy.signal.filter.highpass(filtcopy[ii].data,freq=5000,df=48000,corners=5, zerophase=True) for ii in range(0,len(sthead))]


















