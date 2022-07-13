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
from scipy.signal import sosfiltfilt, butter, detrend, sosfreqz


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
        files = lst[i:i + sliceLen]
        # file1 , file2 , file3  = [*files] #unpack the 3 input files
        # read the sampling frequencies and unpack
        fs1,fs2,fs3 = \
            [*[TdmsFile.read_metadata(f).properties['SamplingFrequency[Hz]'] for f in files]]
        if not fs1 == fs2 == fs3:
            print(f'Sampling Frequencies of input files are different: Files skipped')
            
            continue
        elif fs1 == fs2 == fs3:
            yield lst[i:i + sliceLen]
            
        
        
        
        
def iDAS_timeAvg(files):
    """
    function to process 30 sec tdms DAS files for low frequency strain rate
    
    check to make sure Fs is the same for all 3 files
    concatenate the files together
    taper the combined data
    low pass the data w a cut freq of 1/30 Hz (30 sec)
    pick the middle point?
    
    Parameters
    ----------
    files : tdms DAS data files
    tdms (National Instrument proprietary file type) files for processing

    Returns
    -------
    filt_dat : pandas dataframe
    dataframe of tapered, low-passed data from selected DAS channels
    """
    # create a single channel or range if you don't want to read all channels
    chan = np.r_[10:20]
    
    file1 , file2 , file3  = [*files] #unpack the 3 input files
    # read the sampling frequencies and unpack
    fs1,fs2,fs3 = \
        [*[TdmsFile.read_metadata(f).properties['SamplingFrequency[Hz]'] for f in files]] 
    # if not fs1 == fs2 == fs3:
    #     print(f'Sampling Frequencies of input files are different: Files skipped')
    #     pass
        

    if fs1 == fs2 == fs3 == float(1000):
        # use the open() method so only the file metadata is read initially
        # the returned Tdms object should be used as a context manager to keep
        # the file open and allow specific channel data to be read on demand
        with TdmsFile.open(file1) as f1:
            # use tdms file 
            channel1 = f1["insert group_name"]["insert channel_name"]
            all_channel_data1 = channel1[:]
            data_subset1 = channel1[chan]
            
         with TdmsFile.open(file2) as f2:
             # use tdms file 
             channel2 = f2["insert group_name"]["insert channel_name"]
             all_channel_data2 = channel2[:]
             data_subset2 = channel2[chan]
             
         with TdmsFile.open(file3) as f3:
             # use tdms file 
             channel3 = f3["insert group_name"]["insert channel_name"]
             all_channel_data3 = channel3[:]
             data_subset3 = channel3[chan]
       
        # concatenate the 3 files 
        three_files = np.concatenate((data_subset1,data_subset2,data_subset3), axis = 0)
        # delete large variables and garbage collect
        del data_subset1, data_subset2, data_subset3
        del all_channel_data1, all_channel_data2, all_channel_data3
        gc.collect()
        # detrend the mean and any linear trends from the data chunk 
        detrend_dat = detrend(three_files,axis=0,type='constant')
        detrend_dat = detrend(detrend_dat,axis=0,type='linear')
        # build the window
        cos_win = cosine(len(three_files))
        employ the window/taper
        win_dat = three_files.multiply(cos_win,axis=0)    
        del three_files
        gc.collect()
      
        # build the filter
        Wn = 1/30 # set the low pass frequency ie 30 s here
        sos = butter(2,Wn,output='sos',fs=fs1)
        # view the filter response 
        # w,h = sosfreqz(sos,worN=512,fs=fs1)
        # db = 20*np.log10(np.maximum(np.abs(h), 1e-5))
        # plot it if you want
        # apply the filter
        filt_dat = sosfiltfilt(sos,win_dat)    
        mid_dat = filt_dat[math.trunc(len(filt_dat)/2)]
        # mid_dat = detrend_dat[math.trunc(len(detrend_dat)/2)] # use this if you only want to look at raw data detrended
        return mid_dat

    elif fs1 == fs2 == fs3 == float(10000):
        # use the open() method so only the file metadata is read initially
        # the returned Tdms object should be used as a context manager to keep
        # the file open and allow specific channel data to be read on demand
        with TdmsFile.open(file1) as f1:
            # use tdms file 
            channel1 = f1["insert group_name"]["insert channel_name"]
            all_channel_data1 = channel1[:]
            data_subset1 = channel1[chan]
            
         with TdmsFile.open(file2) as f2:
             # use tdms file 
             channel2 = f2["insert group_name"]["insert channel_name"]
             all_channel_data2 = channel2[:]
             data_subset2 = channel2[chan]
             
         with TdmsFile.open(file3) as f3:
             # use tdms file 
             channel3 = f3["insert group_name"]["insert channel_name"]
             all_channel_data3 = channel3[:]
             data_subset3 = channel3[chan]
        # concatenate the 3 files 
        three_files = np.concatenate((data_subset1,data_subset2,data_subset3), axis = 0)
        # delete large variables and garbage collect
        del data_subset1, data_subset2, data_subset3
        del all_channel_data1, all_channel_data2, all_channel_data3
        gc.collect()
        # detrend the mean and any linear trends from the data chunk 
        detrend_dat = detrend(three_files,axis=0,type='constant')
        detrend_dat = detrend(detrend_dat,axis=0,type='linear')
        # build the window
        cos_win = cosine(len(three_files))
        # employ the window/taper
        win_dat = three_files.multiply(cos_win,axis=0)    
        del three_files
        gc.collect()
        # build the filter
        Wn = 1/30 # set the low pass frequency ie 30 s here
        sos = butter(2,Wn,output='sos',fs=fs1)
        # apply the filter
        filt_dat = sosfiltfilt(sos,win_dat)
        mid_dat = filt_dat[math.trunc(len(filt_dat)/2)]
        # mid_dat = detrend_dat[math.trunc(len(detrend_dat)/2)] # use this if you only want to look at raw data detrended
        return mid_dat

#%% 
# You have to cd into the directory where the data live bc mp.pool doesn't 
# work via interactive interpreters
start = time.time()
dirpath = 'insert_path_data'
os.chdir(dirpath)
postfix='tdms'
file_list = [f for f in sorted(os.listdir(dirpath)) if (postfix in f)]

#change slice_len and number of processors as neccessary
slice_len = 3
num_proc = 12
# start the pool
pool = mp.Pool(processes = num_proc)

# use list comprehension to feed the files from sliceIterator into the iDAS_timeAvg func
proc = [pool.apply_async(iDAS_timeAvg,args=(files,)) for files in sliceIterator(file_list,slice_len)]    
results = [p.get() for p in proc]
print(time.time()-start)
# if you want convert to dataframe and save
results = pd.DataFrame(results)
results.to_csv(r'/set/to/path/for/processed/data',\
                      header=True,index=None,sep=',',mode='a')

