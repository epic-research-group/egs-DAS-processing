#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 13:26:24 2022

@author: spri902
"""

import matplotlib.pyplot as plt
%matplotlib qt5
import numpy as np
import pandas as pd
import ray_tr_fcns as rtf


# input MARMOD data from Table 4.1 in Shearer text
# this is a generic marine seismic model that assumes linear velocity
# gradients within each of the 4 layers
depth = [0.0, 1.5, 6.0, 6.5, 10]
pVel  = [4.5, 6.8, 7.0, 8.0, 8.1]
sVel  = [2.4, 3.75, 3.85, 4.6, 4.7]
dns   = [2.0, 2.8, 2.9, 3.1, 3.1]

# set columns for dataframe and build dataframe
colnms = ['depth', 'pVel', 'sVel', 'dns']
marmod = pd.DataFrame(list(zip(depth,pVel,sVel,dns)), columns=colnms)
# create a vector of 100 evenly spaced ray parameters to loop over
p   = np.linspace(0.1236,0.2217,100)
f = [] # output list for distance and time values

for ii in range(len(p)): # loop over each ray parameter

    x = [0] # for each iter set dist and time lists to zero
    t = [0]
    
    for jj in range(len(marmod.depth)-1): # loop over depths in model
    
        utop = 1/pVel[jj]   # set the slowness at the top of each layer
        ubot = 1/pVel[jj+1] # set the slowness at the bot of each layer
        
        # find the layer thickness
        h = np.diff(marmod.depth)[jj]  
        [dx, dt, irtr] = rtf.layerxt(p[ii], h, utop, ubot) # trace the rays
        
        x += dx # add up each iter of dist and time
        t += dt
        if irtr == 2: # if the ray turns in this layer go to the next iter
            continue
        
    f.append([2*x,2*t,p[ii]]) # collect the dist and time from the inner loop

# convert to a dataframe and sort the data from early to late
fdf = pd.DataFrame(f[::-1],columns=['distance','time','rayP'])
#fdf = pd.DataFrame(2*pd.DataFrame(f[::-1],columns=['dist','time']).values,columns=['dist','time'])
# make the reduced time data   
fdf['redVel'] = fdf.time-fdf.distance/8
# make the tau(p) vector
fdf['tau'] = fdf.time - fdf.rayP*fdf.distance

# plot T(X) curves for part a
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(fdf.distance,fdf.redVel)
plt.xlabel('Distance (km)')        
plt.ylabel('Reduced Travel Time (s)') 
plt.title('Travel Time vs Distance w/ Reduction Velocity of 8 km/s')     
ax.annotate(
    'prograde',
    xy=(20, 0.75), xycoords='data',
    xytext=(-50, 30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'retrograde',
    xy=(42, 1.3), xycoords='data',
    xytext=(-25, 25), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'prograde',
    xy=(50, 1.3), xycoords='data',
    xytext=(30, -30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'caustic',
    xy=(21.8, 1.08), xycoords='data',
    xytext=(30, -30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'caustic',
    xy=(74.7, 1.8), xycoords='data',
    xytext=(30, -30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))


# plot X(p)
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(fdf.rayP,fdf.distance)
plt.title('X(p)')
plt.xlabel('Ray paramter (s/km)')
plt.ylabel('Distance (km)')
ax.annotate(
    'prograde',
    xy=(0.125, 90), xycoords='data',
    xytext=(50, 30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'retrograde',
    xy=(0.1355, 37), xycoords='data',
    xytext=(-35, 25), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'prograde',
    xy=(0.15, 40), xycoords='data',
    xytext=(30, -30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'prograde',
    xy=(0.185, 10), xycoords='data',
    xytext=(-30, 30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'caustic',
    xy=(0.1255, 22.3), xycoords='data',
    xytext=(30, -30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'caustic',
    xy=(0.1436, 74.2), xycoords='data',
    xytext=(30, -30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))

# plot Tau(p)
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(fdf.rayP,fdf.tau)
plt.title('Tau(p)')
plt.xlabel('Ray paramter (s/km)')
plt.ylabel('Tau (p) Intercept')
ax.annotate(
    'prograde',
    xy=(0.1257, 1.109), xycoords='data',
    xytext=(15, 10), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'retrograde',
    xy=(0.1408, 0.725), xycoords='data',
    xytext=(15, 10), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
ax.annotate(
    'prograde',
    xy=(0.1479, 0.301), xycoords='data',
    xytext=(15, 10), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))