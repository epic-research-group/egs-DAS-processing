#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 20:40:47 2022

@author: spri902
"""
import os
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# class structtype():
#     pass
# well = structtype()
wells = pd.read_csv('/home/spri902/Collab_metadata/WellPoints.csv')
wells.columns = wells.columns.str.replace(' ','')
wells.x = wells.x * 0.3048
wells.y = wells.y * 0.3048
wells.z = wells.z * 0.3048

well = {}
# well['I']   = wells.loc[wells['HoleID'] == 'E1-I']
# well['P']   = wells.loc[wells['HoleID'] == 'E1-P']
well['OT']  = wells.loc[wells['HoleID'] == 'E1-OT']
well['OB']  = wells.loc[wells['HoleID'] == 'E1-OB']
well['PST'] = wells.loc[wells['HoleID'] == 'E1-PST']
well['PSB'] = wells.loc[wells['HoleID'] == 'E1-PSB']
well['PDT'] = wells.loc[wells['HoleID'] == 'E1-PDT']
well['PDB'] = wells.loc[wells['HoleID'] == 'E1-PDB']

wn     = ['OT','OB','PST','PSB','PDB','PDT']
ch_bot = [396, 576, 837, 982, 1150, 1329]
max_z  = [59.39, 58.07, 41.01, 59.31, 58.71, 58.17]
ch_in  = ch_bot - np.round(max_z).astype(int)
ch_out = ch_bot + np.round(max_z).astype(int)

chans = np.empty((1728,4),dtype=object)
chans[:] = np.NaN
# chans = [[0 for x in range(4)] for y in range(1728)]
for i in range(6):
    x = well[wn[i]]['Depth(m)'].values
    xp = [well[wn[i]].x.values, well[wn[i]].y.values, well[wn[i]].z.values]
    fp = np.r_[ch_bot[i]:ch_out[i]+1] - ch_bot[i]
    gp = np.r_[ch_in[i]:ch_bot[i]+1] - ch_in[i]
    f = interp1d(x,xp)
    chans[ch_bot[i]:ch_out[i]+1,0] = wn[i]
    chans[ch_in[i]:ch_bot[i]+1,0] = wn[i]
    chans[ch_bot[i]:ch_out[i]+1,1:] = f(fp).T
    chans[ch_in[i]:ch_bot[i]+1,1:] = f(gp).T
    
chancols = ['HoleID','x','y','z']   
chans = pd.DataFrame(chans)
chans.columns = chancols

chans.to_csv(r'/home/spri902/Collab_metadata/DASchannels.txt',\
                       columns=chancols,sep=',')