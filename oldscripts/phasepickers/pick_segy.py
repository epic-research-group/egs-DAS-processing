#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 17:51:20 2022

@author: spri902
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:24:24 2021

@author: spri902
"""

#%%
import pickle
import numpy as np
import pandas as pd
%matplotlib
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from mpl_toolkits.mplot3d import Axes3D
import os
from obspy.io.segy.segy import _read_segy
from obspy import read
import scaleogram as scg
from phasepickers import *
from scipy.signal import *
import seaborn as sns
#from scipy import stats
#from scipy.optimize import curve_fit
#from tsmoothie.smoother import *

geomcols=np.arange(8)
wellcols=np.arange(12)
CASSMGeom=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='Geode',header=0,usecols=geomcols)
CASSMsrc=pd.read_excel('/home/spri902/Collab_metadata/SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx',\
                        sheet_name='CASSM_Channels',header=0,usecols=geomcols)
WellGeom=pd.read_excel('/home/spri902/Collab_metadata/Well_Points.xlsx',header=0,usecols=wellcols)

freq=5       #frequency content (kHz)
unc=0.1    #data uncertainty (ms)

src=[11]
recChan=slice(12,23)
otgeom=[]
#cc_vecsamps=np.array([])
# pstgeom=[]
dfcols=['count','sourceID','receiverID','picktime','std','filename','timestamp']
OT_header=pd.DataFrame([])
OTcorr_picks=pd.DataFrame([])
# PSTcorr_picks=pd.DataFrame([])

OTtrigger_index=np.array([0])  
OTtrigger_ptnl_index=[]
OTtrigger_times=np.array([0])

outlr=[]

# low = 3500
# high = 8000
cnt=1
idx=-1
ncnt=np.empty(shape=(1,2))
ncCols=['count','picktime']
#fl noise ratio filter window length before and after potential picks used to calculate standard deviation       
#N control threshold level to determine if remove the pick by comparing std or rms on both sides of each potential pick  
#%%
# for path,subdir,files in sorted(os.walk('/home/spri902/EGS_Collab/joint_inv/')):
#     for filename in sorted(files):
#         if filename.endswith(".segy"):
    path='/shared/EGS_CASSM/20180426001344'
    for filename in os.listdir(path):
        if filename.endswith(".segy"):
            nsigma = [3,4.25] 
            head = _read_segy(os.path.join(path,filename), headonly=True)
            sthead = read(os.path.join(path,filename))
            times = sthead[0].times()
            srcCh = CASSM_chanDecode(t=times, enc=sthead[93].data)
            srcID = CASSMsrc[CASSMsrc["CytekCh"] == srcCh]
            dt = sthead[0].stats.delta
            fs = np.floor(sthead[0].stats.sampling_rate)
            nyq = np.floor(0.5*(1/dt))
            st = sthead[0].stats.starttime
            et = sthead[0].stats.endtime
            sthead.trim(st, st + (dt*999))
            data = np.stack(t.data for t in sthead.traces).T
            del sthead
            vm = np.percentile(data[:,12:23], 99)
            print("The 99th percentile is {:.0f}; the max amplitude is {:.0f}".format(vm, data.max()))
            plt.figure(figsize=(18,6))
            plt.imshow(data[:,12:23], cmap="RdBu", vmin=-vm, vmax=vm, aspect='auto')
            plt.colorbar()
            plt.show()
            fhdr=head.traces[0].header
            OThy = sthead[recChan]#+sthead[72:75]
            # PSTtr = sthead[51:58]+sthead[59]#+sthead[72:75]
            
            #OThy = OThy.filter('bandpass', freqmin=low,freqmax=high,
                               #corners=4, zerophase=True)
            # PSTtr = PSTtr.filter('bandpass', freqmin=low,freqmax=high,
            #                    corners=4, zerophase=True)            
            OThydata = np.stack(tr.data for tr in OThy.traces).T
            # PSTdata=np.stack(tr.data for tr in PSTtr.traces).T
            #sthead.detrend("linear")
            #sthead.detrend("demean")

            # time of moving avg window (trace length in seconds / some value )
            t_ma = (sthead[0].stats.npts)*(sthead[0].stats.delta)/150
            # number of samples in the moving avg window
            npts_Tma = int(round(t_ma/dt, 0))
            slen = sthead[0].stats.npts  # signal length
            t = np.arange(0, slen / (1/dt), dt)  # signal time vector
            if srcCh in src:
                if srcCh==src[0]:
                    for ii in range(len(OThy)):
                        #trigger the picks
                        OTthresh = threshold(OThy[ii], t_ma, nsigma[0])
                        OTtrigger_index = np.where(aic_parker(
                            OThy[ii][npts_Tma:slen]) > OTthresh[npts_Tma:slen])[0]
                        OTtrigger_index = OTtrigger_index + np.array(npts_Tma)

                        OTtrigger_ptnl_index.append(OTtrigger_index[0])
                        OTrecID = int(
                            float(OThy[ii].stats.seg2['RECEIVER_LOCATION']))
                        tmp=[cnt,srcID.iloc[0,0], OTrecID, t[OTtrigger_index[0]],unc,filename,path[-14:]]
                        otgeom.append([cnt,srcID.iloc[0,0], OTrecID, t[OTtrigger_index[0]],unc,filename,path[-14:]])
                        cnt = cnt+1
                elif srcCh==src[1]:
                    for ii in range(len(OThy)):
                        #trigger the picks
                        OTthresh = threshold(OThy[ii], t_ma, nsigma[1])
                        OTtrigger_index = np.where(aic_parker(OThy[ii][npts_Tma:slen]) > OTthresh[npts_Tma:slen])[0]
                        OTtrigger_index = OTtrigger_index + np.array(npts_Tma)

                        OTtrigger_ptnl_index.append(OTtrigger_index[0])
                        OTrecID = int(
                            float(OThy[ii].stats.seg2['RECEIVER_LOCATION']))
                        tmp=[cnt,srcID.iloc[0,0], OTrecID, t[OTtrigger_index[0]],unc,filename,path[-14:]]
                        otgeom.append([cnt,srcID.iloc[0,0], OTrecID, t[OTtrigger_index[0]],unc,filename,path[-14:]])
                        cnt = cnt+1
            else:
                continue
            OTtrigger_times = np.array(t[OTtrigger_ptnl_index])
            OTtrigger_times = OTtrigger_times.reshape(-1, 1)


            # now run AIC picktimes through the iterative cross correlation picker
                        
            #pwin = [int(x / 2) for x in OTtrigger_ptnl_index]            
            #win = int(np.min(OTtrigger_ptnl_index)/2)
            #ot_vecsamps, lags, _, sample_set=iterative_cc(OThydata, t, OTtrigger_times, win, dt, 300, 300, 25)
            #pwin = int(np.min(OTtrigger_ptnl_index2)/2)
            #pt_vecsamps, lags, _, sample_set=iterative_cc(OThydata, t, OTtrigger_times2, pwin, dt, 300, 300, 25)
            plt.figure()
            seismic_wiggle(OThydata, dt=dt, picktimes=OTtrigger_times,scale=1, color='red', trnormalize=True)
            plt.title(f'Date folder {path.split("/")[-1]} file: {filename} source: {srcCh}')
            # out=plt.ginput(n=-1,timeout=0,show_clicks=True)
            # nout=np.array([round(out[i][0]) for i,j in enumerate(out)])
            # nout=np.array([i+12 for i in nout])
            # nout=np.array([i-(idx*-11) for i in nout])
            # aout=np.array([item[1] for item in out])
            # #ncnt.append(np.array([i-(idx*-11) for i in nout]),np.array([item[1] for item in out]))
            # nout=np.c_[nout,aout]
            # ncnt=np.r_[ncnt,nout]
            #seismic_wiggle(OThydata,dt=dt,picktimes=ot_vecsamps,scale=1,color='blue',trnormalize=True)
            
                      
            OTtrigger_index = np.array([0])
            OTtrigger_ptnl_index = []
            OTtrigger_times = np.array([0])

            idx=idx+1
# sub_picks = pd.DataFrame(ncnt,columns=ncCols)
OTcorr_picks = pd.DataFrame(otgeom, columns=dfcols)
#%%
ndf = OTcorr_picks.merge(sub_picks,on=['count'],how='left',suffixes=('_',''))
ndf['picktime']=ndf['picktime'].fillna(ndf['picktime_']).astype(float)
ndf = ndf.drop('picktime_',axis=1)
OTcorr_picks['picktime'] = ndf['picktime']
#%%
src6 = OTcorr_picks[OTcorr_picks.sourceID==6]
src11= OTcorr_picks[OTcorr_picks.sourceID==11]
plt.figure()
ax1=sns.boxplot(x='receiverID',y='picktime',data=src6)
sns.swarmplot(x='receiverID',y='picktime',data=src6,color="0.25")

#outlrs6 = [y for stat in boxplot_stats(src6['picktime']) for y in stat['fliers']]
plt.title('Source 6')
plt.figure()
ax2=sns.boxplot(x='receiverID',y='picktime',data=src11)
sns.swarmplot(x='receiverID',y='picktime',data=src11,color="0.25")
plt.title('Source 11')



ol6 = find_outlierIQR(src6)
gl6 = find_outlierHampel(src6)
ol11 = find_outlierIQR(src11)
gl11 = find_outlierHampel(src11)
#%% Build E4D Survey File & Save pick data
# OTcorr_picks.to_csv(r'/home/spri902/EGS_Collab/picks/OT_PST_picks.txt',\
#                      header=None,index=None,sep=' ',mode='a')
# sub_picks.to_csv(r'/home/spri902/EGS_Collab/picks/sub_picks.txt',\
#                       header=None,index=None,sep=' ',mode='a')
# CASSMGeom.iloc[12:23].to_csv(r'/home/spri902/EGS_Collab/picks/RECgeometry.txt',\
#                  index=None,sep=' ',mode='a')
# CASSMsrc.to_csv(r'/home/spri902/EGS_Collab/picks/SRCgeometry.txt',\
#                   index=None,sep=' ',mode='a')   
# headerlist=[len(src)+len(OThy), 1]
# headerlist.append(CASSMGeom.iloc[recChan,0].values,CASSMGeom.iloc[recChan,5:].values)
# OT_header = pd.DataFrame(len(src)+len(OThy))
#%% Check straight ray velocity beetween source/receiver
sl_velocity(CASSMsrc.iloc[11:12,5:].values,CASSMGeom.iloc[12:23,5:].values,OTcorr_picks.picktime[0:11].values)
#%% Do some plotting
fig=plt.figure()
for i in np.arange(12,23):
    plt.plot(OTcorr_picks[(OTcorr_picks.receiverID==i) & (OTcorr_picks.sourceID==6)].picktime.values,'.')
fig.legend()
plt.title('Source 6')
fig1=plt.figure()
for i in np.arange(12,23):
    plt.plot(OTcorr_picks[(OTcorr_picks.receiverID==i) & (OTcorr_picks.sourceID==11)].picktime.values,'.')
fig1.legend()
plt.title('Source 11')



fig=plt.figure()
ax = fig.add_subplot(projection='3d')
s1=ax.scatter(WellGeom.iloc[:,3:9:6]/3.28084,WellGeom.iloc[:,4:10:6]/3.28084,WellGeom.iloc[:,5:11:6]/3.28084,\
           label='Well Points',marker='.',color='black',s=2)
s2=ax.scatter(CASSMGeom.iloc[12:23,5],CASSMGeom.iloc[12:23,6],CASSMGeom.iloc[12:23,7],label='OT hydro',marker='o',color='blue')
s3=ax.scatter(CASSMGeom.iloc[69:78:3,5],CASSMGeom.iloc[69:78:3,6],CASSMGeom.iloc[69:78:3,7],label='OT accels',marker='^',color='green')
s4=ax.scatter(CASSMGeom.iloc[60:69:3,5],CASSMGeom.iloc[60:69:3,6],CASSMGeom.iloc[60:69:3,7],label='OB accels',marker='^',color='green')
s5=ax.scatter(CASSMGeom.iloc[51:60:3,5],CASSMGeom.iloc[51:60:3,6],CASSMGeom.iloc[51:60:3,7],label='PST accels',marker='^',color='green')
s6=ax.scatter(CASSMsrc.iloc[0:6,5],CASSMsrc.iloc[0:6,6],CASSMsrc.iloc[0:6,7],label='OB Sources',marker='*',color='red')
s7=ax.scatter(CASSMsrc.iloc[6:12,5],CASSMsrc.iloc[6:12,6],CASSMsrc.iloc[6:12,7],label='PST Sources',marker='*',color='red')
s8=ax.scatter(CASSMsrc.iloc[12:17,5],CASSMsrc.iloc[12:17,6],CASSMsrc.iloc[12:17,7],label='PSB Sources',marker='*',color='magenta')
s9=ax.scatter(CASSMsrc.iloc[17,5],CASSMsrc.iloc[17,6],CASSMsrc.iloc[17,7],label='PDT Source',marker='*',color='red')
s10=ax.scatter(CASSMGeom.iloc[:12,5],CASSMGeom.iloc[:12,6],CASSMGeom.iloc[:12,7],label='PDB hydro',marker='o',color='blue')

ax.set_xlabel('Easting')
ax.set_ylabel('Northing')
ax.set_zlabel('Elevation')
ax.legend()
plt.show()



#%%
#seismic_wiggle(OThydata,dt=dt,picktimes=ot_vecsamps,scale=1,color='blue',trnormalize=True)
#seismic_wiggle(OThydata,dt=dt,picktimes=pt_vecsamps,scale=1,color='green',trnormalize=True)
#scg.cws(t,OThydata[:,10],scales=np.logspace(0,8,num=256,endpoint=True,base=2),wavelet='cmor1.5-1.0', \
#figsize=(14, 7), cmap="viridis", cbar=None, ylabel="Period [seconds]", \
#xlabel="Time [seconds]",yscale='log')
#f,tx,Sxx = spectrogram(sthead[77].data,1/dt)
#plt.pcolormesh(tx,f,Sxx,shading='gouraud')
#plt.ylabel('Frequency [Hz]')
#plt.xlabel('Time [sec]')
#plt.show()



































































