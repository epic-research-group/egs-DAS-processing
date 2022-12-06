#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: spri902
%
% USAGE: err = computeErrorFunction( u1, u0, nSample, lag, norm )
%
% INPUT:
%   u1      = trace that we want to warp; size = (nsamp,1)
%   u0      = reference trace to compare with: size = (nsamp,1)
%   nSample = numer of points to compare in the traces
%   lag     = maximum lag in sample number to search
%   norm    = 'L2' or 'L1' (default is 'L2')
% OUTPUT:
%    err = the 2D error function; size = (nsamp,2*lag+1)
%
% The error function is equation 1 in Hale, 2013. You could umcomment the
% L1 norm and comment the L2 norm if you want on Line 29
%
% Original by Di Yang
% modified by Dylan Mikesell (25 Feb. 2015)
% ported to Python Parker Sprinkle (31, Jan. 2022)
"""
import cupy as cp
import pickle

def computeErrFunc(u1,u0,nSample,lag=80,norm='L2'):
    if lag >= nSample:
        raise ValueError
        print('Lag value must be less than the number of samples in trace')
        
    # Allocate space for error function variable
    err = cp.zeros((nSample,2*lag + 1))
    
    # initial error function calculation
    for ll in cp.arange(-lag,lag + 1): # loop over the lags PS edited for 0 index
        thisLag = ll + lag #P.S. version -> python indexed from zero
        
        for ii in cp.arange(nSample): # loop over the samples
        # skip corners for now we will come back to these
            if (ii + ll >= 0) and (ii + ll <= nSample - 1): #P.S. edited for 0 indexing
                
                diff = u1[ii] - u0[ii + ll] # take sample difference
                
                if norm == "L2":
                    
                    err[ii,thisLag] = diff**2
                    
                elif norm == "L1":
                    
                    err[ii,thisLag] = abs(diff)
    
    # now fix the corners with constant extrapolation
    
    for ll in cp.arange(-lag,lag+1): # loop over lags P.S. edited for 0 index
        thisLag = ll + lag #P.S. edited for 0 index
        
        for ii in cp.arange(nSample): #loop over samples
        
            
            if (ii + ll < 1): # lower left corner (negative lag, early time)
                
                err[ii,thisLag] = err[-ll + 1,thisLag]
            elif (ii + ll > nSample): # upper right corner (pos lag, late time)
                
                err[ii,thisLag] = err[nSample - ll, thisLag]
                
    return err           
    
def accumulateErrFunc(direction,err,nSample,lag=80,b=2):
    """    
    Parameters
    ----------
    direction : integer
        accumulation direction ( dir > 0 = forward in time, dir <= 0 = backward in time).
    err : numpy array from computeErrorFunction output
        the 2D error function; size = (nsamp,2*lag+1).
    nSample : integer
        numer of points to compare in the traces.
    lag : integer, optional
        maximum lag in sample number to search < nSample. The default is 80.
    b : integer, optional
        strain limit (integer value >= 1). The default is 1.

    Returns
    -------
    d = the 2D distance function; size = (nsamp,2*lag+1)
    
    The function is equation 6 in Hale, 2013.
    Original by Di Yang
    Modified by Dylan Mikesell (25 Feb. 2015)
    Ported to Python by Parker Sprinkle (31, Jan. 2022)

    """                      
    nLag = (2*lag) + 1 # num of lags from [-lag:+lag]
    # allocate space for the distance array
    d = cp.zeros((nSample,nLag))
    #--------------------------------------------------------------------------
    # Setup indices based on forward or backward accumulation direction
    #--------------------------------------------------------------------------
    if direction > 0: 
        iBegin = 0 #P.S. edits here for 0 indexing
        iEnd   = nSample-1 #P.S. edits here for 0 indexing
        iInc   = 1
        loop   = cp.arange(iBegin,iEnd + 1,iInc) # P.S. version
    else:
        iBegin = nSample-1 #P.S. edits here for 0 indexing
        iEnd   = 0 #P.S. edits here for 0 indexing
        iInc   = -1
        loop   = cp.arange(iBegin,iEnd - 1,iInc) # P.S. version
    #--------------------------------------------------------------------------
    # Loop through all times in forward or backward direction
    for ii in loop:
        
        # min/max to account for the edges/boundaries
        ji = max(1,min(nSample-1,(ii+1)-iInc)) # i-1 index
        jb = max(1,min(nSample-1,(ii+1)-(iInc * b))) # i-b index
        
        # loop through all lags
        for ll in cp.arange(nLag): # 0 -> 400
            #------------------------------------------------------------------
            #check limits on lag indices
            lMinus1 = ll - 1 # lag at l-1
            
            #if lMinus1 < 1:
            #   lMinus1 = 1
            if lMinus1 < 0: #P.S. version
                lMinus1 = 0
            
            lPlus1 = ll + 1 # lag at l + 1
            
            #if lPlus1 > nLag: original
            if lPlus1 >= nLag: # P.S. version
                lPlus1 = nLag - 1 #D.Y. version
                #lPlus1 = nLag # D.M. version
            #------------------------------------------------------------------
            
            # get distance at lags (ll-1, ll, ll+1)
            distLminus1 = d[jb,lMinus1] # minus:  d( i-b,j-1)
            distL       = d[ji,ll]      # actual: d( i-1,j)
            distLplus1  = d[jb,lPlus1]  # plus:   d( i-b,j+1)
            
            if ji != jb: #equation 10 in Hale (2013)
                #for kb in cp.arange(ji,jb + iInc,-iInc):
                for kb in cp.arange(ji,jb + iInc + 1,-iInc): # P.S. version
                    distLminus1 = distLminus1 + err[kb,lMinus1]
                    distLplus1  = distLplus1  + err[kb,lPlus1]
            
            # Hale equation 6 (if b=1) or 10 (if b>1) after treating boundaries
            d[ii,ll] = err[ii,ll] + min([distLminus1,distL,distLplus1])
    
    return d
    
def backtrackDistFunc(direction, d, err, lmin, b):
    """
    USAGE:    

    Input Parameters
    ----------
    direction : TYPE
        side to start minimization ( dir > 0 = front, dir <= 0 =  back).
    d : TYPE
        the 2D distance function; size = (nsamp,2*lag+1).
    err : TYPE
        the 2D error function; size = (nsamp,2*lag+1).
    lmin : TYPE
        minimum lag to search over.
    b : TYPE
        strain limit (integer value >= 1).

    Returns
    -------
    stbar : vector of integer shifts subject to |u(i)-u(i-1)| <= 1/b.
    
    The function is Equation 2 in Hale ,2013
    Original by Di Yang
    Modified by Dylan Mikesell Dec 2014
    Ported to Python by Parker Sprinkle Jan 2022

    """ 
    nSample = d.shape[0] #num of samples
    nLag    = d.shape[1] # num of lags
    stbar   = cp.zeros((nSample)) # set space for warping function
    
    #--------------------------------------------------------------------------
    # Setup indices based on forward or backward accumulation direction
    #--------------------------------------------------------------------------
    
    if direction > 0:    # FORWARD
        iBegin = 0       # start index P.S. edits here for 0 indexing
        iEnd   = nSample-1 # end index
        iInc   = 1       # increment
    else:                # BACKWARD
        iBegin = nSample-1 # start index
        iEnd   = 0       # stop index
        iInc   = -1      # increment
    
    #--------------------------------------------------------------------------
    # start from the end (front or back)
    ll = cp.argmin(d[iBegin,:]) # min accum dist at front/back based on 'dir'
    
    stbar[iBegin] = ll + lmin # abs value of integer shift
    #--------------------------------------------------------------------------
    # move through all time samples in forward or backward direction
    
    ii = iBegin
    while ii != iEnd:
        
        # min/max for edges/boundaries
        ji = max(1,min(nSample-1,ii + iInc))
        jb = max(1,min(nSample-1,ii + iInc*b))
        
        #----------------------------------------------------------------------
        # check limits on lag indices
        lMinus1 = ll - 1
        
        if lMinus1 < 0: # P.S. edits for 0 indexing
            lMinus1 = 0
        # if lMinus1 < 1: # check lag index is greater than 1
        #     lMinus1 = 1 # make lag = first lag
        
        lPlus1 = ll + 1 # lag at l+1
        
        # if lPlus1 > nLag: # check lag index less than max lag
        #     lPlus1 = nLag # D. M. and D. Y. version
            #lPlus1 = nLag-1
        if lPlus1 >= nLag: # P.S. edited for 0 indexing
            lPlus1 = nLag - 1
        #----------------------------------------------------------------------
        
        # get distance at lags (ll-1, ll, ll+1)
        distLminus1 = d[jb,lMinus1] # minus:  d[ i-b, j-1]
        distL       = d[ji,ll]      # actual: d[ i-1, j]
        distLplus1  = d[jb,lPlus1]  # plus:   d[ i-b,j+1]
        
        if ji != jb: # equation 10 in Hale 2013
        
            ###### what should the end be below (+ 1) or (- 1)#####
            for kb in cp.arange(ji,jb-iInc + 1,iInc): # sum errors over i-1:i-b+1
                
                distLminus1 = distLminus1 + err[kb,lMinus1]
                distLplus1  = distLplus1  + err[kb,lPlus1]
        
        dl = min([distLminus1, distL, distLplus1]) # update min dist to previous sample
        
        if dl != distL: # then ll != ll and we check forward and backward
            if dl == distLminus1:
                ll = lMinus1
            else: # dl == lPlus1
                ll = lPlus1
        
        # assume ii = ii - 1
        ii = ii + iInc # previous time sample
        stbar[ii] = ll + lmin  # absolute integer of lag
        
        # now move to correct time index , if smoothing difference over many
        # time samples using 'b'
        if ll == lMinus1 or ll == lPlus1: # check edges to see about b values
            if ji != jb: # if b > 1 then need to move more steps
                for kb in cp.arange(ji,(jb - iInc) + 1,iInc):
                    ii = ii + iInc
                    stbar[ii] = ll + lmin 
    return stbar

def computeDTWerr(Aerr, u, mLag):
    """
    Parameters
    ----------
    Aerr : TYPE
        error MATRIX (eq 13 in Hale 2013)
    u : TYPE.
       warping function (samples) VECTOR 
    mLag : TYPE
        value of maximum lag (samples) SCALER

    Returns
    -------
    error
    Written by Dylan Mikesell Feb 2015
    Ported to Python by Parker Sprinkle Jan 2022

    """
    npts = len(u) 
    
    if cp.shape(Aerr,1) != npts: #### check this #### 
        
        print('Error matrix dimensions are incorrect: Check inputs')
        Aerr = Aerr.T
    
    error = 0 #initialize error
    
    # accumulate error 
    for ii in cp.arange(1,npts+1):
        idx = mLag + 1 + u[ii] #index of lag
        error = error + Aerr[ii,idx]
    
    return error
    

def save_file(obj,filename):
    try:
        with open(f'{filename}.pickle',"wb") as f:
            pickle.dump(obj,f,protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error during pickling object",ex)
def load_object(filename):
    try:
        with open(filename,"rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during pickling object",ex)
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    