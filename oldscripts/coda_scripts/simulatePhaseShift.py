#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 16:26:27 2022

@author: spri902
"""

import numpy as np
import matplotlib
%matplotlib
import matplotlib.pyplot as plt

npts = 2000
A = 0.7
f = 5
w = 2.0*np.pi*f
dt = 1/500
t = np.linspace(0,dt*npts,npts)
phi = np.pi/6
dphi = np.pi/2

u0 = 1 + A*np.cos(w*t) + 
plt.plot(t,u0)
u1 = A*np.cos(w*t - dphi)
#%% 
