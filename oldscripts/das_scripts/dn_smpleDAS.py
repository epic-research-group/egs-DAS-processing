#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 22:36:36 2022

@author: spri902
"""

from nptdms import TdmsFile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib
from scipy.signal import *
import obspy


fs=10000
dt=1/fs
tdms_file = TdmsFile.read("/shared/EGS_iDAS/20180521/Sigma_Exp1_180521000032.tdms")
for group in tdms_file.groups():
    group_name = group.name
    for channel in group.channels():
        channel_name = channel.name
        # Access dictionary of properties:
        properties = channel.properties
        # Access numpy array of data for channel:
        data = channel[:]
        t=np.arange(0,dt*len(data),dt)
        data = data-np.mean(data)
        # Access a subset of data
        #data_subset = channel[100:200]
        ff = butter(4, 1,analog=False,output='sos',fs=fs) #default is lowpass
        ss = sosfiltfilt(ff, data)
        plt.plot(ss)
