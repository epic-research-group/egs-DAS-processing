#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 21:52:53 2020

@author: spri902
"""

import pickle
import pandas as pd
import numpy as np

c2sourceWpickdf = pickle.load( open( "/Users/spri902/OneDrive - PNNL/Desktop/RT_picker/results/subTer_picks/sourceWpickdf", "rb" ) )
c2sourceSpickdf = pickle.load( open( "/Users/spri902/OneDrive - PNNL/Desktop/RT_picker/results/subTer_picks/sourceSpickdf", "rb" ) )
subterGeom=pd.ExcelFile('/Users/spri902/Desktop/SubTER_BCD/SubTER/subTER_Source_rec.xlsx')
sourceW=pd.read_excel(subterGeom,'SourceW',header=1)
sourceS=pd.read_excel(subterGeom,'SourceS',header=1)
rcoor=pd.read_excel(subterGeom,'PhonesN',header=1)
# e4dhead=pd.read_excel(subterGeom,'Header')

ndata=[]

ffname='SubTer'# name of file
freq=8       #frequency content (kHz)
unc=0.125    #data uncertainty (ms)

# Write E4D text survey file
# number of source-receiver pairs 1
# source-reciever pair id, x y z, frequency content (kHZ)
# blank line
# number of data, source id, receiver id, ttime (ms), uncertainty (ms)

# get source and receiver coordinates
frames=[sourceW,sourceS]
scoor=pd.concat(frames,join='outer')
scoor=np.array(scoor)
rcoor=np.array(rcoor)
ch=rcoor[:,1].astype(int)
ns=len(scoor[:,1])
nr=len(rcoor[:,1])

ep=np.zeros_like(scoor[:,1]).T
for ii in range(ns):
    ind=np.all((scoor[:,1]==scoor[ii,1], scoor[:,2]==scoor[ii,2],scoor[:,3]==scoor[ii,3]),axis=0)
    ep[ind]=ii
rp=np.zeros_like(rcoor[:,1]).T
for ii in range(nr):
    ind=np.all((rcoor[:,2]==rcoor[ii,2], rcoor[:,3]==rcoor[ii,3],rcoor[:,4]==rcoor[ii,4]),axis=0)
    rp[ind]=ii
rp=rp+np.max(ep+1)
sr=nr+ns

# c2sourceWpickdf=c2sourceWpickdf.explode('picktime').reset_index(drop=True)
# c2sourceWpickdf=c2sourceWpickdf.loc[c2sourceWpickdf['picktime'] != 0.]
# c2sourceSpickdf=c2sourceSpickdf.explode('picktime').reset_index(drop=True)
# c2sourceSpickdf=c2sourceSpickdf.loc[c2sourceSpickdf['picktime'] != 0.]
# ndata=len(c2sourceSpickdf.index)+len(c2sourceWpickdf.index)
ndata=len(c2sourceSpickdf.explode('picktime'))+len(c2sourceWpickdf.explode('picktime'))
pickdf=pd.concat([c2sourceWpickdf,c2sourceSpickdf])
# time1=c2sourceWpickdf['picktime']
# time2=c2sourceSpickdf['picktime']
# ttime=pd.concat([time1,time2]).reset_index()


##
with open('subTer_survey.txt', 'a') as fid:
    fid.write('%-4d 1\n' % sr)
    
    # sources
    for ii in range(ns):
    
        fid.write('%-4d %-7.3f %-7.3f %-6.2f %-3d\n' % (ii,scoor[ii,1],scoor[ii,2],scoor[ii,3],freq))
    # receivers
    for jj in range(nr):
        fid.write('%-4d %-7.3f %-7.3f %-6.2f %-3d\n' % (jj+ns,rcoor[jj,2],rcoor[jj,3],rcoor[jj,4],freq))
    # travel times

    fid.write('\n')
    fid.write('%-8d\n' % ndata)
    cnt=0
    
    for i in range(ns):
        for j in range(nr):
            fid.write('%-8d %-4d %-4d %-7.7f %-7.3f\n' % (cnt+1,pickdf.iloc[i,1][0],pickdf.iloc[i,2]['count'][j],pickdf.iloc[i,3][j][0],unc))
            cnt=cnt+1
            
            
            
            
            
            
            
            
            
            