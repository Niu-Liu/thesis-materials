#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 21:51:56 2016

@author: Neo
"""

import numpy as np

def ICRF_IERS():
    datd = '../list/'
    datf = ['icrf2-defining.dat', 'icrf2-non-vcs.dat', 'icrf2-vcs-only.dat']
    row = [20, 23, 20]
    
    IERS = np.array([], dtype=str)
    ICRF = np.array([], dtype=str)
    for i in range(len(datf)):
        icrf, iers = np.loadtxt(datd+datf[i], dtype=str, skiprows=row[i],\
                                delimiter='  ', usecols=(0,1), unpack=True) 
        IERS, ICRF = np.hstack((IERS, iers)), np.hstack((ICRF, icrf))
        
    return IERS, ICRF
    
#IERS, ICRF = ICRF_IERS()
#print IERS.size