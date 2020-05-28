#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 12:47:50 2016

@author: Neo

Generate machine-readable versions of our tables
"""

import numpy as np
import time
from SourceDesignation import ICRF_IERS

## ICRF and IERS designation for sources.
IERS, ICRF = ICRF_IERS()

## time series statistics
sta_fil = '../results/sou.sta'
Souo = np.loadtxt(sta_fil, dtype=str, usecols=(0,))
sesn = np.loadtxt(sta_fil, dtype=int, usecols=(1,))
obsp = np.loadtxt(sta_fil, usecols=(4,))

## variance data
var_fil = '../results/variance.dat'
Souv = np.loadtxt(var_fil, dtype=str, usecols=(0,))
Stda, WASa, Stdd, WASd = \
    np.loadtxt(var_fil, usecols=(3,5,7,9), unpack=True)
    
#ICRF2 defining sources data file
icrf2f = '../data/icrf2-defining.dat'
defn = np.loadtxt(icrf2f, dtype=str, skiprows=20, usecols=(2,))

## rank lists
ran_fil = ['../results/OV.rank', '../results/GR.rank']
## two final list
ind = np.array([323, 294])

## data format suitable for latex
tex_fil = '../results/table.tex'
ftex = open(tex_fil, 'w')
print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
            time.localtime(time.time())), file=ftex)

for i in range(ind.size):
    num = ind[i]
    out_fil = '../results/Sou%d'%num + '.dat'
    fout = open(out_fil, 'w')
    
    print('##Sou%d'%num, file=ftex)
    
    Sour = np.loadtxt(ran_fil[i], dtype=str, usecols=(0,)) 
    V, V_E, V_RA, V_Dec = \
        np.loadtxt(ran_fil[i], usecols=(1,2,4,7), unpack=True) 
    
## only subset are needed.
    Sous, Vs, V_Es, V_RAs, V_Decs =\
        Sour[:num], V[:num], V_E[:num], V_RA[:num], V_Dec[:num]
## normalized.
    NVs = Vs/V_Es    
    
    for j in range(num):
        k = np.where(Souv==Sous[j])[0][0]
        l = np.where(Souo==Sous[j])[0][0]
        m = np.where(IERS==Sous[j])[0][0]
        
#        if Sous[j] in defn:
#            ty = 'D'
#        else:
#            ty = 'N'
            
        print('%s  %s'%(ICRF[m], Sous[j]) + \
            ('  %6.2f'+'  %6.1f'*2)%(NVs[j], V_RA[j], V_Dec[j]) +\
            ('  %6.2f'*4)%(Stda[k], Stdd[k], WASa[k], WASd[k]) +\
            ('  %6.1f  %4d')%(obsp[l], sesn[l]), file=fout) 
#            +\
#            '  %s'%ty
## data format used in future for CDS            
#   1-  8  I8    ---            Source    TERS Designations
#  11- 16  F6.6  ---            NLD       Normalized linear drift
#  19- 21  D3    0.001mas/yr    ldRA      Linear drift of RA.cos(DEC.)
#  24- 26  D3    0.001mas/yr    e_ldRA    Formal uncertainty
#  29- 31  D3    0.001mas/yr    ldDEC     Linear drift of DEC.
#  34- 36  D3    0.001mas/yr    e_ldDEC   Formal uncertainty
#  39- 43  F5.3  mas            WStdRA    Weighted Standard deviation of RA.cos(DEC.)
#  46- 50  F5.3  mas            WStdDEC   Weighted Standard deviation of DEC.
#  53- 57  F5.5  mas            WAStdRA   Weighted Allan Standard deviation of RA.cos(DEC.)
#  60- 64  F5.5  mas            WStdDEC   Weighted Allan Standard deviation of DEC.        
#  67- 70  F4.1  yr             T_obs     Observation span
#  73- 77  D4    ---            N_ses     Number of sessions
#  80      I1    ---            type      Source type in ICRF2 catalogue     

        if j < 10:
            print('%s  &%s'%(ICRF[m], Sous[j]) + \
            ('  &%6.2f'+'  &%6.1f'*2)%(NVs[j], V_RA[j], V_Dec[j]) +\
            ('  &'+'  &%6.2f'*2+'  &'+'  &%6.2f'*2)%(Stda[k], Stdd[k], WASa[k], WASd[k]) +\
            ('  &%6.2f  &%4d\\\\')%(obsp[l], sesn[l]), file=ftex) 
#            +\
#            '  %s \\\\'%ty