#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 07:43:02 2022

@author: spri902
"""
from itertools import groupby
from nptdms import TdmsFile
import pandas as pd
import numpy as np
import time
import os
#%%
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
        
def all_equal(iterable):
    g = groupby(iterable)
    return next(g,True) and not next(g,False)
 #%%       
start = time.time()
dirpath = '/data1/parker/EGS_iDAS/20180526/'
os.chdir(dirpath)
postfix='tdms'
file_list = [f for f in sorted(os.listdir(dirpath)) if (postfix in f)]
slice_len = 3

skipFiles = []
for files in sliceIterator(file_list,slice_len):
    file1 , file2 , file3  = [*files] #unpack the 3 input files
    # read the sampling frequencies and unpack
    fs1,fs2,fs3 = \
        [*[TdmsFile.read_metadata(f).properties['SamplingFrequency[Hz]'] for f in files]] 
    
    # print(fs1,fs2,fs3)
    # print(file1,file2,file3)
    if not fs1 == fs2 == fs3:
        # print(f'Sampling Frequencies of input files are different: {files} skipped')
        skipFiles.append(files)
print(time.time()-start)
skipFiles = pd.DataFrame(skipFiles)
skipFiles.to_csv(r'/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/lpFilter/skipFlies.txt',\
                      header=True,index=None,sep=',',mode='a')