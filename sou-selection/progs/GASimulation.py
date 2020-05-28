#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:34:11 2016

@author: Neo

Simulation of Galactic Abberation.

"""

import numpy as np
cos = np.cos
sin = np.sin
import time
from RotationFit import MeanWeiRotation

def simulation(alp, det):
## Notice !!!
## unit for RA and DE are rad.    
    
##    For GC
##    alp0 = radians(266.25)
##    det0 = radians(-29)
## For the North Pole
    alp0 = np.deg2rad(0)
    det0 = np.deg2rad(90)
    A = 5
    
    d1 = A*cos(alp0)*cos(det0)
    d2 = A*sin(alp0)*cos(det0)
    d3 = A*sin(det0)

## a simulated proper motion field.   
    ua = -d1*sin(alp)          + d2*cos(alp)
    ud = -d1*cos(alp)*sin(det) - d2*sin(alp)*sin(det) + d3*cos(det)
             
    return ua, ud

## some path variables
res_dir = '../results/'
dat_fil = ['OV.rank', 'GR.rank']
#dat_fil = ['OV1.rank', 'GR1.rank']
res = ['OVsim.dat', 'GRsim.dat']

for d in range(len(dat_fil)):
    fout = open(res_dir + res[d],'w')
    print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
                               time.localtime(time.time())), file=fout)
    print('## Data file: %s'%dat_fil[d], file=fout)
    print('''## Global Spin Vector of the Simulated proper motion field.
##data format
##Num   w_tol   w_x   w_y   w_z
##            (Expectation | Std)
##unit: uas/yr
############################################################################''', file=fout)

    Sou = np.loadtxt(res_dir + dat_fil[d], dtype=str, usecols=(0,)) 
    V, V_err, RA, V_RA, Va_E, Dec, V_Dec, Vd_E = \
        np.loadtxt(res_dir + dat_fil[d], usecols=list(range(1,9)), unpack=True) 
## Unit: deg -> rad
    RA = np.deg2rad(RA)
    Dec = np.deg2rad(Dec)
    
    N= len(Sou)
    
    for i in range(100,len(Sou)+1):
        sig = np.ones(i)
        ra, dec = RA[:i], Dec[:i]
        ua, ud = simulation(ra, dec)
        
        w, sig, corrcoef = MeanWeiRotation(ua, ud, ra, dec)
        wx, wy, wz = w
        sigwx, sigwy, sigwz = sig
        wt = np.sqrt(sum(w**2))
        sigw = np.sqrt(sum(sig**2))
    
        print(("%3d"+"    %8.4f"*8)%(i, wt, sigw, wx, sigwx, wy, sigwy, wz, sigwz), file=fout)
        
fout.close()
print('Done!')