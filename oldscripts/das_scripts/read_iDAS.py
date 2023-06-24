#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 10:56:05 2022

@author: spri902
"""
import os
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
from scipy.integrate import cumulative_trapezoid
from boreholes import parse_surf_boreholes
#%%
surf_wells = ['OT', 'OB', 'PSB', 'PST', 'PDB', 'PDT']

attr_map = {'OT': ['otDepths', 'otTemps'], 'OB': ['obDepths', 'obTemps'],
            'PDB': ['pdbDepths', 'pdbTemps'], 'PDT': ['pdtDepths', 'pdtTemps'],
            'PSB': ['psbDepths', 'psbTemps'], 'PST': ['pstDepths', 'pstTemps']}
well_dict = parse_surf_boreholes('/home/spri902/Collab_metadata/WellPoints.csv')

#%%
def iDAS_timeAvg(filename):
    
    tdms_file = TdmsFile(filename).as_dataframe()
    
    dat=tdms_file.mean().values # mean of data matrix (time avg down columns)
    # dat=tdms_file.apply(abs).mean().values
    # dat=tdms_file.mad().values # mean abs deviation of matrix
    # dat = abs(tdms_file - tdms_file.median())
    # dat = dat.median().values # median abs deviation of matrix
#     #print(filename)
#     #print(time.time() - start)
    
    return dat
# You have to cd into the directory where the data live bc mp.pool doesn't 
# work via interactive interpreters

start = time.time()
dirpath = '/data1/parker/EGS_iDAS/20180525/'
os.chdir(dirpath)
postfix='tdms'
file_list = [f for f in os.listdir(dirpath) if (postfix in f)]
num_proc = 8
pool = mp.Pool(processes = num_proc)
proc = [pool.apply_async(iDAS_timeAvg,args=[filename]) for filename in sorted(os.listdir(dirpath))]
results = [p.get() for p in proc]
dat = pd.DataFrame(list(map(np.ravel, results))).transpose()
# fd = [name.split("_") for name in file_list]
# fl = [fd[file][2].split(".") for file in range(len(fd))]
# fl = [el[0] for el in fl]
# dates = [datetime.strptime(d,'%y%m%d%H%M%S') for d in sorted(fl)]
# xlims = mdates.date2num(dates)
dat.to_csv(r'/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/timeAvgMean_data/20180525.txt',\
                      header=True,index=None,sep=',',mode='a')
print(time.time()-start)
#%% Put separate dat files back together
os.chdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/timeAvgMean_data/')
datFiles = sorted(os.listdir('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/timeAvgMean_data'))
df_all = pd.concat((pd.read_csv(f,header=0,sep=',') for f in datFiles if f.endswith('.txt')),axis=1)
# ndf = [df_all[col].str.split(' ',expand=True) for col in df_all.columns]
# ndf = pd.concat(ndf,axis=1)
df_all.to_pickle('maystim22_26_combined')
# and to read the pickle
df_all = pd.read_pickle('/home/spri902/EGS_Collab/4850/results/maystim/processed_DAS/timeAvgMean_data/maystim22_26_combined')
#mask the OT well data and remove column means from the remaining data
OTchans  = np.r_[337:456]
OBchans  = np.r_[518:635]
PSTchans = np.r_[796:879]
PSBchans = np.r_[923:1042]
PDBchans = np.r_[1091:1210]
PDTchans = np.r_[1271:1388]
 
# df_all = df_all - df_all.mean()
# Do integration in time and space to get back to displacement from strain rate
# integrate down columns for strain rate to strain
# df_all = pd.DataFrame(cumulative_trapezoid(df_all.values,axis=0,initial=0))
# df_all = pd.DataFrame(cumulative_trapezoid(df_all.values,axis=1,initial=0))

  


df_all.mask((df_all.iloc[0:OTchans[-1]+1,:].notna()),inplace=True)
# df_all.mask((df_all.iloc[0:OTchans[0],:].notna()),inplace=True) 
df_all.mask((df_all.iloc[OTchans[-1]+1:OBchans[0],:].notna()),inplace=True)  
df_all.mask((df_all.iloc[OBchans[-1]+1:PSTchans[0],:].notna()),inplace=True) 
df_all.mask((df_all.iloc[PSTchans[-1]+1:PSBchans[0],:].notna()),inplace=True)  
df_all.mask((df_all.iloc[PSBchans[-1]+1:PDBchans[0],:].notna()),inplace=True) 
df_all.mask((df_all.iloc[PDBchans[-1]+1:PDTchans[0],:].notna()),inplace=True)  
df_all.mask((df_all.iloc[PDTchans[-1]+1:1728].notna()),inplace=True)
#%%
wn     = ['OT','OB','PST','PSB','PDB','PDT']
nfile_list = sorted(os.walk('/data1/parker/EGS_iDAS'))
nfile_list = nfile_list[1:6]
# file_list = file_list[1:]
nfile_list = [group[2] for group in nfile_list]
nfile_list = [item for sublist in nfile_list for item in sublist]
# [file_list.append(f) for f in nfile_list]
fd = [name.split("_") for name in nfile_list]
fl = [fd[file][2].split(".") for file in range(len(fd))]
fl = [el[0] for el in fl]
dates = [datetime.strptime(d,'%y%m%d%H%M%S') for d in sorted(fl)]
xlims = mdates.date2num(dates)  
chans=np.linspace(0,df_all.shape[0] - 1,df_all.shape[0]).astype(int)
ch_bot = [396, 576, 837, 982, 1150, 1329]
#%%
# with PdfPages(f'EGS_Collab_iDAS_May242018_stimulation_scaled.pdf') as pdf: 
         
fig1,ax = plt.subplots()
ax2 = ax.secondary_yaxis("right")
vm = np.nanpercentile(df_all,50)
if vm < 0:
    img=ax.pcolormesh(dates,chans,df_all,cmap='RdBu',vmin=2*vm,vmax=vm+abs(vm))
if vm > 0:
    img=ax.pcolormesh(dates,chans,df_all,cmap='RdBu',vmin=vm/2,vmax=vm*2)
plt.gca().invert_yaxis()
#ax.hlines(min(xlims),max(xlims),color='k')
#ax.text(xlims[50],400,"OT well", color='black', ha="center", va="center", bbox =dict(facecolor='none', boxstyle='square, pad=0.3', lw=2))
ax.xaxis_date()
ax2.set_yticks(ch_bot)
#ax2.set_ylabels(wn)
plt.colorbar(img,orientation="horizontal")
#plt.gcf().autofmt_xdate()
plt.title(f'iDAS data {dates[0]} to {dates[-1]}')
#     pdf.savefig(fig1)
#     plt.close()
#%% Stack the DAS data in each well

OBdown = OBchans[0:math.trunc(len(OBchans)/2)]
OBup   = np.flipud(OBchans[math.trunc(len(OBchans)/2)+1:])
OBdat  = (df_all.iloc[OBdown,:].values + df_all.iloc[OBup,:].values)/2

PSBdown = PSBchans[0:math.trunc(len(PSBchans)/2)]
PSBup   = np.flipud(PSBchans[math.trunc(len(PSBchans)/2)+1:])
PSBdat  = (df_all.iloc[PSBdown,:].values + df_all.iloc[PSBup,:].values)/2

PSTdown = PSTchans[0:math.trunc(len(PSTchans)/2)]
PSTup   = np.flipud(PSTchans[math.trunc(len(PSTchans)/2)+1:])
PSTdat  = (df_all.iloc[PSTdown,:].values + df_all.iloc[PSTup,:].values)/2

PDBdown = PDBchans[0:math.trunc(len(PDBchans)/2)]
PDBup   = np.flipud(PDBchans[math.trunc(len(PDBchans)/2)+1:])
PDBdat  = (df_all.iloc[PDBdown,:].values + df_all.iloc[PDBup,:].values)/2

PDTdown = PDTchans[0:math.trunc(len(PDTchans)/2)]
PDTup   = np.flipud(PDTchans[math.trunc(len(PDTchans)/2)+1:])
PDTdat  = (df_all.iloc[PDTdown,:].values + df_all.iloc[PDTup,:].values)/2

ch_List = [OBdown,PSBdown,PSTdown,PDBdown,PDTdown]
df_List = [OBdat,PSBdat,PSTdat,PDBdat,PDTdat]
w_List = wn[1:]

for i in range(len(ch_bot)-1):
    fig1,ax = plt.subplots()
    # ax2 = ax.secondary_yaxis("right")
    vm=np.nanpercentile(df_List[i],50)
    if vm < 0:
        img=ax.pcolormesh(dates,ch_List[i],df_List[i],cmap='RdBu',vmin=2*vm,vmax=vm+abs(vm))
    if vm > 0:
        img=ax.pcolormesh(dates,ch_List[i],df_List[i],cmap='RdBu',vmin=vm/2,vmax=vm*2)
    plt.gca().invert_yaxis()
    #ax.hlines(min(xlims),max(xlims),color='k')
    #ax.text(xlims[50],400,"OT well", color='black', ha="center", va="center", bbox =dict(facecolor='none', boxstyle='square, pad=0.3', lw=2))
    ax.xaxis_date()
    ax2.set_yticks(ch_bot)
    #ax2.set_ylabels(wn)
    plt.colorbar(img,orientation="horizontal")
    #plt.gcf().autofmt_xdate()
    plt.title(f'iDAS_Well_{w_List[i]}_{dates[0]} to {dates[-1]}')
    pickle.dump(ax, open(f'iDAS_Well_{w_List[i]}','wb'))
    plt.close()

for i in range(len(ch_bot)-1):
    fig = pickle.load(open(f'iDAS_Well_{w_List[i]}','rb'))