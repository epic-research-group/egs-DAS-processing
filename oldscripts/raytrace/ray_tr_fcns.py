#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 13:19:35 2022
Functions for ray tracing through the MARMOD oceanic crust model

functions below from Shearer Introduction to Seismology
@author: spri902
"""

# LAYERXT calculates dx and dt for a ray in a layer with a linear
#   velocity gradient.  This is a highly modified version of a
#   Fortran subroutine in Chris Chapman's WKBJ program.
#
# Inputs:   p     =  horizontal slowness
#           h     =  layer thickness
#           utop  =  slowness at top of layer
#           ubot  =  slowness at bottom of layer
# Returns:  dx    =  range offset
#           dt    =  travel time
#           irtr  =  return code
#                 = -1,  zero thickness layer
#                 =  0,  ray turned above layer
#                 =  1,  ray passed through layer
#                 =  2,  ray turned in layer, 1 leg counted in dx,dt
import numpy as np
def layerxt(p, h, utop, ubot):
  # returns [dx, dt, irtr]
  
  if p >= utop:
     return [0.0, 0.0, 0]
  elif h==0:
     return [0.0, 0.0, -1]
  else:
    b = 1.0/(ubot*h) - 1.0/(utop*h)
    def eta(u): return np.sqrt(u**2 - p**2)
    def x(u): return eta(u)/(u*b*p)
    def t(u): return (np.log((u + eta(u))/p))/b
    if utop == ubot:
       return [h*p/eta(utop), h*utop**2/eta(utop), 1 ]
    elif p>=ubot:
       return [x(utop), t(utop), 2]
    else:
       return [x(utop) - x(ubot), t(utop) - t(ubot), 1]



# RTCOEF calculates P/SV reflection/transmission coefficients
# for interface between two solid layers, based on equations 5.40
# (p. 144-145) of Aki and Richards (2nd edition).

#  Inputs:    vp1     =  P-wave velocity of layer 1 (top layer)
#  (real)     vs1     =  S-wave velocity of layer 1
#             den1    =  density of layer 1
#             vp2     =  P-wave velocity of layer 2 (bottom layer)
#             vs2     =  S-wave velocity of layer 2
#             den2    =  density of layer 2
#             hslow   =  horizontal slowness (ray parameter)
#  Returns:   rt(0)   =  down P to P up     (refl)
#  (complex)  rt(1)   =  down P to S up     (refl)
#             rt(2)   =  down P to P down   (tran)
#             rt(3)   =  down P to S down   (tran)
#             rt(4)   =  down S to P up     (refl)
#             rt(5)   =  down S to S up     (refl)
#             rt(6)   =  down S to P down   (tran)
#             rt(7)   =  down S to S down   (tran)
#             rt(8)   =    up P to P up     (tran)
#             rt(9)  =    up P to S up     (tran)
#             rt(10)  =    up P to P down   (refl)
#             rt(11)  =    up P to S down   (refl)
#             rt(12)  =    up S to P up     (tran)
#             rt(13)  =    up S to S up     (tran)
#             rt(14)  =    up S to P down   (refl)
#             rt(15)  =    up S to S down   (refl)

# NOTE:  All input variables are real.
#        All output variables are complex!
#        Coefficients are not energy normalized.

def RTCOEF(vp1, vs1, den1, vp2, vs2, den2, hslow):
   
   rt = np.zeros([16,], dtype=complex)

   alpha1 = complex(vp1, 0.)
   beta1 = complex(vs1, 0.)
   rho1 = complex(den1, 0.)
   alpha2 = complex(vp2, 0.)
   beta2 = complex(vs2, 0.)
   rho2 = complex(den2, 0.)
   p = complex(hslow, 0.)

   cone = complex(1., 0.)
   ctwo = complex(2., 0.)

   si1 = alpha1*p
   si2 = alpha2*p
   sj1 = beta1*p
   sj2 = beta2*p

   ci1 = np.sqrt(cone - si1**2)
   ci2 = np.sqrt(cone - si2**2)
   cj1 = np.sqrt(cone - sj1**2)
   cj2 = np.sqrt(cone - sj2**2)

   term1 = (cone - ctwo*beta2*beta2*p*p)
   term2 = (cone - ctwo*beta1*beta1*p*p)
   a = rho2*term1 - rho1*term2
   b = rho2*term1 + ctwo*rho1*beta1*beta1*p*p
   c = rho1*term2 + ctwo*rho2*beta2*beta2*p*p
   d = ctwo*(rho2*beta2*beta2 - rho1*beta1*beta1)
   E = b*ci1/alpha1 + c*ci2/alpha2
   F = b*cj1/beta1 + c*cj2/beta2
   G = a - d*ci1*cj2/(alpha1*beta2)
   H = a - d*ci2*cj1/(alpha2*beta1)
   DEN = E*F + G*H*p*p

   trm1 = b*ci1/alpha1 - c*ci2/alpha2
   trm2 = a + d*ci1*cj2/(alpha1*beta2)
   rt[0] = (trm1*F - trm2*H*p*p)/DEN              #refl down P to P up

   trm1 = a*b + c*d*ci2*cj2/(alpha2*beta2)
   rt[1] = (-ctwo*ci1*trm1*p)/(beta1*DEN)       #refl down P to S up

   rt[2] = ctwo*rho1*ci1*F/(alpha2*DEN)         #trans down P to P down

   rt[3] = ctwo*rho1*ci1*H*p/(beta2*DEN)        #trans down P to S down

   trm1 = a*b + c*d*ci2*cj2/(alpha2*beta2)
   rt[4] = (-ctwo*cj1*trm1*p)/(alpha1*DEN)      #refl down S to P up

   trm1 = b*cj1/beta1 - c*cj2/beta2
   trm2 = a + d*ci2*cj1/(alpha2*beta1)
   rt[5] = -(trm1*E - trm2*G*p*p)/DEN             #refl down S to S up

   rt[6] = -ctwo*rho1*cj1*G*p/(alpha2*DEN)      #trans down S to P down

   rt[7] = ctwo*rho1*cj1*E/(beta2*DEN)          #trans down S to S down


   trm1 = b*ci1/alpha1 - c*ci2/alpha2
   trm2 = a + d*ci2*cj1/(alpha2*beta1)
   rt[10] = -(trm1*F + trm2*G*p*p)/DEN            #refl up P to P down

   trm1 = a*c + b*d*ci1*cj1/(alpha1*beta1)
   rt[11] = (ctwo*ci2*trm1*p)/(beta2*DEN)       #refl up P to S down

   rt[8] = ctwo*rho2*ci2*F/(alpha1*DEN)         #trans up P to P up

   rt[9] = -ctwo*rho2*ci2*G*p/(beta1*DEN)      #trans up P to S up

   trm1 = a*c + b*d*ci1*cj1/(alpha1*beta1)
   rt[14] = (ctwo*cj2*trm1*p)/(alpha2*DEN)      #refl up S to P down

   trm1 = b*cj1/beta1 - c*cj2/beta2
   trm2 = a + d*ci1*cj2/(alpha1*beta2)
   rt[15] = (trm1*E + trm2*H*p*p)/DEN             #refl up S to S down

   rt[12] = ctwo*rho2*cj2*H*p/(alpha1*DEN)      #trans up S to P up

   rt[13] = ctwo*rho2*cj2*E/(beta1*DEN)         #trans up S to S up

   return rt      



# RTCOEFSH calculates SH reflection/transmission coefficients
# for interface between two solid layers, based on equations 5.33
# (p. 139) of Aki and Richards (2nd edition).

#  Inputs:    vs1     =  S-wave velocity of layer 1 (top layer)
#  (real)     den1    =  density of layer 1
#             vs2     =  S-wave velocity of layer 2
#             den2    =  density of layer 2
#             hslow   =  horizontal slowness (ray parameter)
#  Returns:   rt(0)   =  down S to S up     (refl)
#  (complex)  rt(1)   =  down S to S down   (tran)
#             rt(2)  =   up S to S up       (tran)
#             rt(3)  =   up S to S down     (refl)

# NOTE:  All input variables are real.
#        All output variables are complex!
#        Coefficients are not energy normalized.

def RTCOEFSH(vs1, den1, vs2, den2, hslow):
   
   rt_sh = np.zeros([4,], dtype=complex)

   beta1 = complex(vs1, 0.)
   rho1 = complex(den1, 0.)
   beta2 = complex(vs2, 0.)
   rho2 = complex(den2, 0.)
   p = complex(hslow, 0.)

   cone = complex(1., 0.)

   sj1 = beta1*p
   sj2 = beta2*p
   cj1 = np.sqrt(cone - sj1**2)
   cj2 = np.sqrt(cone - sj2**2)

   D = rho1*beta1*cj1 + rho2*beta2*cj2
   rt_sh[0] = (rho1*beta1*cj1 - rho2*beta2*cj2)/D
   rt_sh[3] = -rt_sh[0]
   rt_sh[1] = 2.*rho1*beta1*cj1/D
   rt_sh[2] = 2.*rho2*beta2*cj2/D

   return rt_sh



# GETAUX returns auxilary fault plane strike, dip & rake,
# given strike,dip,rake of main fault plane.
def GETAUX(strike1, dip1, rake1):
   degrad = 180./3.1415927
   s1 = strike1/degrad
   d1 = dip1/degrad
   r1 = rake1/degrad

   d2 = np.arccos(np.sin(r1)*np.sin(d1))

   sr2 = np.cos(d1)/np.sin(d2)
   cr2 = -np.sin(d1)*np.cos(r1)/np.sin(d2)
   r2 = np.arctan2(sr2, cr2)

   s12 = np.cos(r1)/np.sin(d2)
   c12 = -1./(np.tan(d1)*np.tan(d2))
   s2 = s1 - np.arctan2(s12, c12)

   strike2 = s2*degrad
   dip2 = d2*degrad
   rake2 = r2*degrad

   if (dip2 > 90.):
      strike2 = strike2 + 180.
      dip2 = 180. - dip2
      rake2 = 360. - rake2
   if (strike2 > 360.): strike2 = strike2 - 360.

   return strike2, dip2, rake2

def my_annotate(ax, s, xy_arr=[], *args, **kwargs):
  ans = []
  an = ax.annotate(s, xy_arr[0], *args, **kwargs)
  ans.append(an)
  d = {}
  try:
    d['xycoords'] = kwargs['xycoords']
  except KeyError:
    pass
  try:
    d['arrowprops'] = kwargs['arrowprops']
  except KeyError:
    pass
  for xy in xy_arr[1:]:
    an = ax.annotate(s, xy, alpha=0.0, xytext=(0,0), textcoords=an, **d)
    ans.append(an)
  return ans