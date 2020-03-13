#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 10:27:03 2016

@author: Neo

ICRF offset due to the Galactic aberration.

"""

import numpy as np
cos = np.cos
sin = np.sin
import time
from RotationFit import RotationFit
#from RotationFit import MeanWeiRotation

def simulation(alp, det):
## Notice !!!
## unit for RA and DE are rad.    
    
##    For GC
    alp0 = np.deg2rad(266.25)
    det0 = np.deg2rad(-29)
### For the North Pole
#    alp0 = np.deg2rad(0)
#    det0 = np.deg2rad(90)
    A = 5
    
    d1 = A*cos(alp0)*cos(det0)
    d2 = A*sin(alp0)*cos(det0)
    d3 = A*sin(det0)

## a simulated proper motion field.   
    ua = -d1*sin(alp)          + d2*cos(alp)
    ud = -d1*cos(alp)*sin(det) - d2*sin(alp)*sin(det) + d3*cos(det)
             
    return ua, ud

## output
flog = open('../logs/galactic_aberration.log', 'a')
print>>flog, time.strftime('%Y-%m-%d %H:%M:%S Begins!',time.localtime(time.time()))
## data file.
datf = '../data/qsocom.dat'    
## read data
RAg, DEg = np.loadtxt(datf, usecols=(1,2), unpack=True)    ## unit: "
#RAi, DEi = np.loadtxt(datf, usecols=(5,6), unpack=True)    ## unit: "
#dRA, dDe = RAg -RAi, DEg - DEi
RAg, DEg = np.deg2rad(RAg/3600), np.deg2rad(DEg/3600)      ## unit: rad
Flag     = np.loadtxt(datf, usecols=(13,),  dtype=str )
Frot     = np.loadtxt(datf, usecols=(14,),   dtype=int ) 

pmRA, pmDE = simulation(RAg, DEg)

indd = np.where(Flag=='D')[0]
indn = np.where(Flag=='N')[0]
ind0 = np.where(Frot==0)[0]
ind1 = np.where(Frot==1)[0]
ind3 = np.where(Frot==3)[0]

label = ['All', 'defining', 'Non-defining', 'frame-fixed']
ind = [np.arange(RAg.size), \
       indd, \
       indn  \
       ]
for i in range(3):
    print>>flog, 'subset '+label[i]+' :'
    ra, de = np.take(RAg, ind[i]),  np.take(DEg,  ind[i])
    ua, ud = np.take(pmRA, ind[i]), np.take(pmDE, ind[i])
    eu = np.ones(ua.size)
    
    print>>flog, 'Weighted: Normal.'
    w, sig, corrcoef = RotationFit(ua, ud, eu, eu, ra, de)
    [wx,   wy,  wz] = w
    [ewx, ewy, ewz] = sig

    print>>flog, 'the spin are(uas*yr-1):\n :\n', \
        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
    print>>flog, '   correlation coefficients are:\n', corrcoef
    print>>flog, 'from J2000.0-2015.0, the cumulative orientation are(uas):\n :\n', \
        '    &$%.3f \pm %.3f$'*3%(wx*15, ewx*15, wy*15, ewy*15, wz*15, ewz*15)
    
## fixed frame ones    
print>>flog, 'subset 260 fixed frame: '

ra  = np.hstack((np.take( RAg, ind3), np.take(RAg, ind1)))
de  = np.hstack((np.take( DEg, ind3), np.zeros(ind1.size)))
ua  = np.hstack((np.take(pmRA, ind3), np.take(pmRA, ind1)))
ud  = np.hstack((np.take(pmDE, ind3), np.take(pmDE, ind1)))
eu = np.ones(ua.size)

    
print>>flog, 'Weighted: Normal.'
w, sig, corrcoef = RotationFit(ua, ud, eu, eu, ra, de)
[wx,   wy,  wz] = w
[ewx, ewy, ewz] = sig

print>>flog, 'the orientation/spin are(uas*yr-1):\n :\n', \
        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
print>>flog, '   correlation coefficients are:\n', corrcoef
print>>flog, 'from J2000.0-2015.0, the cumulative orientation are(uas):\n :\n', \
        '    &$%.3f \pm %.3f$'*3%(wx*15, ewx*15, wy*15, ewy*15, wz*15, ewz*15)
     
print 'Done!'