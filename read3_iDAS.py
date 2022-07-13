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
    # PSTchans = np.r_[796:879]
    # PSBchans = np.r_[923:1042]
    # PDBchans = np.r_[1091:1210] # these are the channel nums of the DAS I want
    PDBchan = 1121
    # PDTchans = np.r_[1271:1388]
    
    file1 , file2 , file3  = [*files] #unpack the 3 input files
    # read the sampling frequencies and unpack
    fs1,fs2,fs3 = \
        [*[TdmsFile.read_metadata(f).properties['SamplingFrequency[Hz]'] for f in files]] 
    # if not fs1 == fs2 == fs3:
    #     print(f'Sampling Frequencies of input files are different: Files skipped')
    #     pass
        

    if fs1 == fs2 == fs3 == float(1000):
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
        # cos_win = cosine(len(three_files))
        # employ the window/taper
        # win_dat = three_files.multiply(cos_win,axis=0)    
        del three_files
        gc.collect()
      
        # build the filter
        # Wn = 1/30 # set the low pass frequency ie 30 s here
        # sos = butter(2,Wn,output='sos',fs=fs1)
        # view the filter response 
        # w,h = sosfreqz(sos,worN=512,fs=fs1)
        # db = 20*np.log10(np.maximum(np.abs(h), 1e-5))
        # apply the filter
        # filt_dat = sosfiltfilt(sos,win_dat)    
        # mid_dat = filt_dat[math.trunc(len(filt_dat)/2)]
        mid_dat = detrend_dat[math.trunc(len(detrend_dat)/2)]
        return mid_dat

    elif fs1 == fs2 == fs3 == float(10000):
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
        # cos_win = cosine(len(three_files))
        # employ the window/taper
        # win_dat = three_files.multiply(cos_win,axis=0)    
        del three_files
        gc.collect()
        # build the filter
        # Wn = 1/30 # set the low pass frequency ie 30 s here
        # sos = butter(2,Wn,output='sos',fs=fs1)
        # apply the filter
        # filt_dat = sosfiltfilt(sos,win_dat)
        # mid_dat = filt_dat[math.trunc(len(filt_dat)/2)]
        mid_dat = detrend_dat[math.trunc(len(detrend_dat)/2)]
        return mid_dat

#%% 
# You have to cd into the directory where the data live bc mp.pool doesn't 
# work via interactive interpreters
start = time.time()
dirpath = '/data1/parker/EGS_iDAS/20180526/'
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
print(time.time()-start)
results = pd.DataFrame(results)
results.to_csv(r'/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/20180526_r.txt',\
                      header=True,index=None,sep=',',mode='a')

#%% Put separate dat files back together
os.chdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/')
datFiles = sorted(os.listdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter'))
df_raw = pd.concat((pd.read_csv(f,header=0,sep=',') for f in datFiles if f.endswith('r.txt')),axis=0)
df_all = pd.concat((pd.read_csv(f,header=0,sep=',') for f in datFiles if f.endswith('n.txt')),axis=0)
df_all.reset_index(drop=True,inplace=True)
df_raw.reset_index(drop=True,inplace=True)
# ndf = [df_all[col].str.split(' ',expand=True) for col in df_all.columns]
# ndf = pd.concat(ndf,axis=1)
# df_all.to_pickle('maystim22_26_combined_r')
wn     = ['OT','OB','PST','PSB','PDB','PDT']
nfile_list = sorted(os.walk('/data1/parker/EGS_iDAS'))
nfile_list = nfile_list[1:]
# file_list = file_list[1:]
nfile_list = [group[2] for group in nfile_list]
nfile_list = [item for sublist in nfile_list for item in sublist]
# [file_list.append(f) for f in nfile_list]
fd = [name.split("_") for name in nfile_list]
fl = [fd[file][2].split(".") for file in range(len(fd))]
fl = [el[0] for el in fl]
dates = [datetime.strptime(d,'%y%m%d%H%M%S') for d in sorted(fl)]
dates = dates[1:-1]
xlims = mdates.date2num(dates)  
chans=np.linspace(0,df_all.shape[0] - 1,df_all.shape[0]).astype(int)
ch_bot = [396, 576, 837, 982, 1150, 1329]
# remove_dates = pd.read_csv('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/skipFlies.txt',sep=',',header=None)

ind2rem = [90,91,1571,1572,3670,3671,4274,4275,6298,6299]
for index in sorted(ind2rem,reverse=True):
    del dates[index]
dates = pd.DataFrame(dates)
