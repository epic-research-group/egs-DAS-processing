#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

iterative cross-correlation

data= data matrix (num_samples,num_traces)
times= time vector for the trace/stream
aicpicktimes= estimated arrival times from above AIC routines 
window for searching during cross correlation in SAMPLES
sampling interval =1/fs
acc val= how many samples should the cc be able to shift 
max shift = parameter for clamping lags ( it lets you move the traces within certain limits)
max iterations = max number of times iterate for optimizing the lag shift between traces

"""

import numpy as np

def cross_correlation(dat_1, dat_2):
    
    dat_1m = dat_1 - np.nanmean(dat_1)
    dat_2m = dat_2 - np.nanmean(dat_2)
    
    dat_1f = np.fft.fft(dat_1m)
    dat_2f_conj = np.conj(np.fft.fft(dat_2m))
    norm_factor = np.sqrt(np.sum(dat_1m**2)*np.sum(dat_2m**2))
        
    return np.real(np.fft.fftshift(np.fft.ifft(dat_1f * dat_2f_conj)/norm_factor))
    
def find_cc_lag(cc_dat):
    
    num_samples = cc_dat.size // 2
    
    lag_vec = np.arange(-num_samples, num_samples)
    
    return lag_vec[np.nanargmax(np.abs(cc_dat))]

def hilbert_envelope(dat):
    
    num_samples = dat.size 
    ind = np.fft.fftfreq(num_samples)
    
    unit_step = np.heaviside(ind, 0.5)
       
    return np.abs(np.fft.ifft(2*np.fft.fft(dat)*unit_step))

def pilot_trace(dat, centers, win):
    
    num_samples, num_traces = dat.shape
    
    dat_sum = np.zeros((2*win, 1))
    for idx in range(num_traces):
        dat_sum += dat[centers[idx]-win:centers[idx]+win, idx]
        
    return dat_sum/num_traces

def data_subset(dat, centers, win):
    
    num_samples, num_traces = dat.shape
    
    dat_sub = np.zeros((2*win, num_traces))
    for idx in range(num_traces):
        dat_sub[:, idx] = dat[int(centers[idx]-win):int(centers[idx]+win), idx]

    return dat_sub

def cross_correlation_reps(dat, pilot_dat):
    
    num_samples, num_traces = dat.shape
    lags = np.zeros((num_traces, 1), dtype = np.int16)
    for idx in range(num_traces):
        cc = cross_correlation(dat[:, idx], pilot_dat)
        lags[idx] = find_cc_lag(cc)
        
    return lags
         

def signal_to_noise(dat, centers, win):
    
    num_samples, num_traces = dat.shape
    
    background = np.sqrt(np.nanmean(dat**2))
    sn_vec = np.zeros((num_traces, 1))
    for idx in range(num_traces):
        
        sn_vec[idx] = \
            np.sqrt(np.nanmean(dat[int(centers[idx]-win):int(centers[idx]+win), idx]**2))
    
    return sn_vec/background

#This function ensures the lags stay within a user specified "max_shift" criterion    
def lag_clamping(lags, max_shift): 
    
    ind = np.argwhere(lags > max_shift)
    lags[ind] = np.sign(lags[ind])*max_shift
    return lags
    
    
def iterative_cc(dat, t_vec, init_picks, win, dt, acc_val, max_shift, max_iter):
    # checks the receivers with arrival picks and filters out levels with no pick information (noisy traces)
    available_lvls = np.argwhere(np.isfinite(init_picks.flatten() == 1)).flatten()
    num_avail_lvls = available_lvls.size
    
    dat_subset = dat[:, available_lvls]
    
    # compute the envelope function to remove the opposite polarity effect
    dat_subset_envelope = np.zeros_like(dat_subset)
    for idx in range(num_avail_lvls): 
        dat_subset_envelope[:, idx] = hilbert_envelope(dat_subset[:, idx])
    
    # input picks are in time, this sends to samples
    sample_set = (init_picks[available_lvls] // dt).astype(np.int16) 
    sample_set_original = sample_set.copy()
    
    # compute the SNR of all receiver levels using the picks
    # and takes out the traces with the highest SNR
    sn_ratio = signal_to_noise(dat_subset_envelope, sample_set, 50) 
    # this is the reference level where we assume the initial picks will be 
    # accurate enough for relative adjusting of picks on other levels
    reference_lvl = np.argmax(sn_ratio)
    
    lags = np.zeros_like(sample_set, dtype = np.int16)   
    
    
    iteration = 0
    
    while iteration < max_iter:
        
        # taking a subset of data around the picks
        dat_refresh = data_subset(dat_subset_envelope, sample_set, win)  
        # pilot trace is the stack of waveforms from all receiver levels
        pilot = np.nanmean(dat_refresh, axis=1)
        
        # do the cross correlation of all waveforms w the pilot trace
        lags_local = cross_correlation_reps(dat_refresh, pilot)
        # update the lags
        lags = lags+lags_local        
        lags = lag_clamping(lags, max_shift)
        sample_set = sample_set_original + lags
        
        if np.sum(np.abs(lags_local)) < (num_avail_lvls // 2):
            
            break # if sum of the lags_local is < half the num of receivers, stop
        
        iteration += 1        
    
    #bringing all receiver levels to the initial pick of the reference level
    reference_lvl_diff = \
        np.abs(sample_set[reference_lvl] - sample_set_original[reference_lvl])
    
    if reference_lvl_diff > acc_val:
        lags -= lags[reference_lvl]
        sample_set = sample_set_original + sample_set[reference_lvl_diff]
        
    ind = np.argwhere(np.abs(lags) > acc_val)
    lags[ind] = 0
    sample_set[ind] = sample_set_original[ind]
    
    #return t_vec[sample_set], lags, np.sum(np.abs(lags_local)), iteration
    return t_vec[sample_set], lags, np.sum(np.abs(lags_local)), sample_set