##!/usr/bin/env python3
## -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:22:14 2017

@author: Neo

Divide the QSOs into 2 sets:  northern, southern declinations.

"""

import numpy as np
import time
from RotationFit3 import RotationFit, RotationFitB, RotationFitC
## function
def Orientation_calc(Condition, Flog):
    ra = np.extract(Condition, RA)
    de = np.extract(Condition, DE)
    d_ra = np.extract(Condition, D_RA)
    d_de = np.extract(Condition, D_DE)
    err_ra = np.extract(Condition, ERR_RA)
    err_de = np.extract(Condition, ERR_DE)
    print('Number of the Sample: %d'%ra.size, file = Flog)
##  Orientation.
    w, sig, corrcoef = RotationFit(d_ra, d_de, err_ra, err_de, ra, de)
    [wx,   wy,  wz] = w*1000
    [ewx, ewy, ewz] = sig*1000
    print('Weighted: Normal. Orientation.', file = FLOG)
    print('the orientation/spin are(uas):\n :\n', \
       '    &$%+3.0f \pm %3.0f$'*3%(wx, ewx, wy, ewy, wz, ewz),\
       file = Flog)
    print('   correlation coefficients are:\n', corrcoef, file = Flog)
## consider a possible bias in declination.    
    print('Weighted: Normal. Consider a bias in declination.', file = Flog)
    w, sig, corrcoef = RotationFitB(d_ra, d_de, err_ra, err_de, ra, de)
    [wx,   wy,  wz, b]  = w*1000
    [ewx, ewy, ewz, eb] = sig*1000
    print('the orientation/spin are(uas):\n :\n', \
       '    &$%+3.0f \pm %3.0f$'*3%(wx, ewx, wy, ewy, wz, ewz),\
       file = Flog)
    print('the bias in declination(mas):\n :\n', \
       '    &$%+3.0f \pm %3.0f$'%(b, eb), file = Flog)
    print('   correlation coefficients are:\n', corrcoef, file = Flog)
## main body
FLOG = open('../logs/qsocom.log', 'a')
print(time.strftime('%Y-%m-%d %H:%M:%S Begins!',\
    time.localtime(time.time())), file = FLOG)
## Load data
RA, DE, D_RA, D_DE, ERR_RA, ERR_DE =\
    np.loadtxt('../data/qsocom.dat',\
        usecols = (1, 2, 9, 10, 11, 12), unpack = True)
RA = np.deg2rad(RA/3600)       ## Unit: rad
DE = np.deg2rad(DE/3600)       ## Unit: rad
## Divide data into 2 sets: Nouthern and Southern declinations.
## For Nouthern hemisphere
con = DE > 0   ## Nouthern
print('## For Nouthern hemisphere:', file = FLOG)
Orientation_calc(con, FLOG)
## For Southern hemisphere
con = DE < 0   ## Southern
print('## For Southern hemisphere:', file = FLOG)
Orientation_calc(con, FLOG)
con = DE != 0  ## All sources
print('## For All sources:', file = FLOG)
Orientation_calc(con, FLOG)