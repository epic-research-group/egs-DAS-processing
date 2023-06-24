#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 21:35:10 2022

@author: spri902
"""

import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import convolve

def ricker(f, length=0.2, dt=1/48000):
    t = np.arange(-length/2, (length-dt)/2, dt)
    y = (1.0 - 2.0*(np.pi**2)*(f**2)*(t**2)) * np.exp(-(np.pi**2)*(f**2)*(t**2))
    return t, y
 
f = 3000 # A low wavelength of 25 Hz
t, w = ricker(f)

rho1 = 2900
v1 = 6000
rho2 = 3100
v2 = 7000

Rt = np.ones(len(t))*(rho2*v2 - rho1*v1) / (rho2*v2 + rho1*v1)

Tr = convolve(Rt,w,mode='same',method='fft')