#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 10:03:01 2022

@author: spri902
"""

import os
import gc
import math
from nptdms import TdmsFile
import pandas as pd
import numpy as np
import multiprocessing as mp
import time
import matplotlib.dates as mdates
from datetime import datetime as datetime
import matplotlib
%matplotlib
import matplotlib.pyplot as plt 
import pickle
from scipy.signal.windows import cosine
from scipy.signal import sosfiltfilt, butter, detrend
from scipy.integrate import cumulative_trapezoid
from boreholes import parse_surf_boreholes

#%% Functions for reading and processing borehole DAS data from EGS Collab Experiment 1

def sliceIterator(lst,sliceLen):
    """
    # function to generate a sliding window of slices for extracting "sliceLen" 
    files at a time from the main directory
    Parameters
    ----------
    lst : list
    lst = list of time sorted files
    sliceLen : int
    number of files to slice from the list

    Yields
    ------
    sliding window selection of DAS files for input into processing function    
    """
    for i in range(len(lst) - sliceLen + 1):
        yield lst[i:i + sliceLen]
        
def iDAS_timeAvg(files):
    """
    Parameters
    ----------
    files : tdms DAS data files
    tdms (National Instrument proprietary file type) files for processing

    Returns
    -------
    filt_dat : pandas dataframe
    dataframe of tapered, low-passed data from selected DAS channels
    """
    # PSTchans = np.r_[796:879]
    # PSBchans = np.r_[923:1042]
    # PDBchans = np.r_[1091:1210] # these are the channel nums of the DAS I want
    PDBchan = 1120
    # PDTchans = np.r_[1271:1388]
    
    file1 , file2 , file3  = [*files] #unpack the 3 input files

    # read them as dataframes with the built-in pandas method in nptdms
    tdms_file1 = TdmsFile(file1).as_dataframe().iloc[:,PDBchan]
    tdms_file2 = TdmsFile(file2).as_dataframe().iloc[:,PDBchan]
    tdms_file3 = TdmsFile(file3).as_dataframe().iloc[:,PDBchan]
    # concatenate the 3 files 
    three_files = pd.concat((tdms_file1,tdms_file2,tdms_file3), axis = 0)
    # delete large variables and garbage collect
    del tdms_file1,tdms_file2,tdms_file3
    gc.collect()
    # detrend the mean and any linear trends from the data chunk 
    detrend_dat = detrend(three_files,axis=0,type='constant')
    detrend_dat = detrend(detrend_dat,axis=0,type='linear')
    # build the window
    cos_win = cosine(len(three_files))
    # employ the window/taper
    win_dat = three_files.multiply(cos_win,axis=0)    
    # build the filter
    Wn = 1/30 # set the critical frequency ie 30 s here
    sos = butter(2,Wn,output='sos')
    # apply the filter
    filt_dat = sosfiltfilt(sos,win_dat)

    
    return filt_dat
#%% 
# You have to cd into the directory where the data live bc mp.pool doesn't 
# work via interactive interpreters
start = time.time()
dirpath = '/data1/parker/EGS_iDAS/20180522/'
os.chdir(dirpath)
postfix='tdms'
file_list = [f for f in sorted(os.listdir(dirpath)) if (postfix in f)]
slice_len = 3
num_proc = 12
pool = mp.Pool(processes = num_proc)

# for files in sliceIterator(file_list,slice_len):
#     proc = pool.apply_async(iDAS_timeAvg,args=[files])
proc = [pool.apply_async(iDAS_timeAvg,args=(files,)) for files in sliceIterator(file_list,slice_len)]    
results = [p.get() for p in proc]
filt_dat = pd.DataFrame(list(map(np.ravel, results))).transpose()
# fd = [name.split("_") for name in file_list]
# fl = [fd[file][2].split(".") for file in range(len(fd))]
# fl = [el[0] for el in fl]
# dates = [datetime.strptime(d,'%y%m%d%H%M%S') for d in sorted(fl)]
# xlims = mdates.date2num(dates)
filt_dat.to_csv(r'/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/timeAvgMean_data/20180525.txt',\
                      header=True,index=None,sep=',',mode='a')
print(time.time()-start)

