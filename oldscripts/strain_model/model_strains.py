#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 04:43:28 2022

@author: spri902
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Analytical solutions found in Zehnder
x1 = np.linspace(0,3,401)
x2 = np.linspace(0,1.5,201)
X1,X2 = np.meshgrid(x1,x2)
a = 1
Z = X1 + 1j*X2

s11 = np.real(Z / np.sqrt(Z**2 - a**2))\
    - X2*np.imag(1/np.sqrt(Z**2 - a**2)\
                 - Z**2/(Z**2-a**2)**(3/2)) - 1

s22 = np.real(Z / np.sqrt(Z**2 - a**2))\
    + X2*np.imag(1/np.sqrt(Z**2 - a**2)\
                 - Z**2/(Z**2-a**2)**(3/2))

E = 1e10
nu = 0.25
e11 = (1+nu)/E *(s11 - nu*(s11+s22))

plt.subplots(1,3,figsize=(10,3))

ax=plt.subplot(1,3,1)
p1=plt.pcolor(x1,x2,s11,vmin=-1,vmax=1,cmap='seismic')
ax.set_aspect('auto')
plt.colorbar(p1)

plt.subplot(1,3,2)
p2=plt.pcolor(x1,x2,s22,vmin=-1,vmax=1,cmap='seismic')
plt.colorbar(p2)

plt.subplot(1,3,3)
p3=plt.pcolor(x1,x2,e11,vmin=-1e-10,vmax=1e-10,cmap='seismic')
plt.colorbar(p3)

plt.tight_layout()
plt.show()

plt.subplots()
plt.plot(x1,e11[25,:])
plt.ylabel('Strain (exx)')
plt.xlabel('Distance along (x / a)')