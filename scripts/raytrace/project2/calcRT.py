#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 14:48:19 2022
Programing project 2
ESS 512 Intro to Seismology
Winter 2022

This script calculates and plots reflection / transmission coefficients for the P/SV 
system for an interface between 2 solid layers based on equations from Aki and 
Richards text. 

The first plot is an R/T amplitude vs ray angle plot for 3 of the 16 coefficients 
generated from the function RTCOEFF that correspond to reflected S and P phases 
and a transmitted P phase. 

The second plot is an R/T phase vs ray angle plot for 3 of the 16 coefficients 
generated from the function RTCOEFF that correspond to reflected S and P phases 
and a transmitted P phase.

@author: spri902
"""
import numpy as np
import pandas as pd
import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import ray_tr_fcns as rtf

# Read in the PREM model data
colnames = ['Depth','Radius','Vp','Vs','density','Qmu','Qk','pressure']
prem = pd.read_csv('~/scripts/raytrace/prem.txt',delim_whitespace=True,\
                   header=0,names=colnames)
# Get wavespeeds and density above and below mantle - core boundary
Vp_m       = prem['Vp'][45]
Vs_m       = prem['Vs'][45]
rho_m      = prem['density'][45]
Vp_c       = prem['Vp'][46]
Vs_c       = 0.000001
rho_c      = prem['density'][46]
# Build ray angle vector
ray_angles = np.linspace(0,90,91)
ray_angles = np.deg2rad(ray_angles)
uS         = 1/Vs_m

# set variable for function output
output = []
# Loop over ray angles 
for angle in ray_angles:
    
    p  = uS*np.sin(angle) # horizontal slowness for this angle
    
    rt = rtf.RTCOEF(Vp_m, Vs_m, rho_m, Vp_c, Vs_c, rho_c, p)
    # get the down S and P up refl
    rone = rt[4]
    # get the down S and S up refl
    rtwo = rt[5]
    # get the down S and P trans
    tone = rt[6]

    #collect the output for each of the 3 R/T phases
    output.append([rone,rtwo,tone])
# Take the real part for the amplitudes and the imaginary part for the phase
amps  = np.abs(output)
phase1 = np.degrees(np.angle([p1[0] for p1 in output]))
phase2 = np.degrees(np.angle([p1[1] for p1 in output]))
phase3 = np.degrees(np.angle([p1[2] for p1 in output]))
# Do some plotting 
# amp vs ray angle 
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(np.rad2deg(ray_angles),amps[:,0],'k',label = 'Reflected P wave')
ax.plot(np.rad2deg(ray_angles),amps[:,1],'r',label = 'Reflected S wave')
ax.plot(np.rad2deg(ray_angles),amps[:,2],'g',label = 'Transmitted P wave')
plt.xlabel('Ray Angle (deg)')        
plt.ylabel('R/T Coefficient Amplitude') 
plt.title('RTCOEFF Amplitudes vs Ray Angle')  
plt.legend()
# phase vs ray angle
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(np.rad2deg(ray_angles),phase1,'k',label = 'Reflected P wave')
ax.plot(np.rad2deg(ray_angles),phase2,'r',label = 'Reflected S wave')
ax.plot(np.rad2deg(ray_angles),phase3,'g',label = 'Transmitted P wave')
plt.xlabel('Ray Angle (deg)')        
plt.ylabel('R/T Coefficient Phase Angle') 
plt.title('RTCOEFF Phase Angle vs Ray Angle')
plt.legend()

 
 
  
