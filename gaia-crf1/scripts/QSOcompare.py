##!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 16:49:07 2016

@author: Neo

Cross Identification between GAIA_qso and icrf2
"""

import numpy as np
import time
from RotationFit import RotationFit, RotationFitB, RotationFitC
import matplotlib.pyplot as plt
#from RotationAndGlide import RGFit
sin = np.sin
cos = np.cos
pi = np.pi

def icrf2def():
    ICRFd = np.loadtxt('../data/icrf2-defining.dat', dtype=str, skiprows=20, \
                        usecols=(1,))
    return ICRFd

def read_icrf2():
    ICRF, sDE = np.loadtxt('../data/icrf2.dat', dtype=str, usecols=(0,7), unpack=True)
    RAh, RAm, RAs, DEd, DEm, DEs, e_RAs, e_DEs, C \
    = np.loadtxt('../data/icrf2.dat', usecols=list(range(4,13,1)), unpack=True)
    DEsig = np.ones(len(sDE))
    for i in range(len(DEsig)):
        if sDE[i][0] == '-':
            DEsig[i] = -1
    RA = 15*((RAh*60+RAm)*60+RAs)             ##   s -> arcsec
    DE = DEsig*((np.abs(DEd)*60+DEm)*60+DEs)  ##   arcsec
    e_RA = e_RAs*15*1000                      ##   s -> mas
    e_DE = e_DEs*1000                         ##   arcsec -> mas
    
    return ICRF, RA, DE, e_RA, e_DE, C
    
def read_qso():
    ICRF2 = np.loadtxt('../data/qso.dat', dtype=str, usecols=(10,), delimiter='|')
    RAdeg, e_RAdeg, DEdeg, e_DEdeg, Cg \
    = np.loadtxt('../data/qso.dat', usecols=list(range(3,8)), unpack=True, delimiter='|')
    Frot = np.loadtxt('../data/qso.dat', usecols=(11,), dtype=int, delimiter='|')
    
    RAg = RAdeg*3600        ##  deg -> arcsec
    e_RAg = e_RAdeg         ##  mas
    DEg = DEdeg*3600        ##  deg -> arcsec
    e_DEg = e_DEdeg         ##  mas
    
    return ICRF2, RAg, DEg, e_RAg, e_DEg, Cg, Frot
    
flog = open('../logs/qsocom.log', 'a')
print(time.strftime('%Y-%m-%d %H:%M:%S Begins!',time.localtime(time.time())), file=flog)
        
ICRFi, RAi, DEi, e_RAi, e_DEi, Ci = read_icrf2()
ICRFg, RAg, DEg, e_RAg, e_DEg, Cg, Frot = read_qso()

#print RAi[0], DEi[0], e_RAi[0], e_DEi[0]
#print RAg[0], DEg[0], e_RAg[0], e_DEg[0]

ind1 = []  ## for gaia
ind2 = []  ## for icrf2
num = 0
N = len(ICRFg)

for j in range(N):
    i = np.where(ICRFi==ICRFg[j])[0]
    if len(i):
        ind2.append(i[0])
        ind1.append(j)
        num += 1
    else:
        print('No corresponding entry for %s are found in icrf2.'%ICRFg[j], file=flog)
            
    print('%d/%d '%(j+1,N)+ ' is completed.')

print(': %d/%d entries are found. Completed!'%(num, N), file=flog)
 
icrf2d = list(icrf2def())
   
ICRFg = np.take(ICRFg, ind1)
RAg = np.take(RAg, ind1)
DEg = np.take(DEg, ind1)
e_RAg = np.take(e_RAg, ind1) 
e_DEg = np.take(e_DEg, ind1)
Cg = np.take(Cg, ind1)
Frot = np.take(Frot, ind1)   

RAi = np.take(RAi, ind2)
DEi = np.take(DEi, ind2)
e_RAi = np.take(e_RAi, ind2)*np.cos(np.deg2rad(DEi/3600))
e_DEi = np.take(e_DEi, ind2)    
Ci = np.take(Ci, ind1)

## offset, gaia - icrf2
d_RA = (RAg - RAi)*np.cos(np.deg2rad(DEg/3600.0))*1.0e3         ##  mas
d_DE = (DEg - DEi)*1.0e3                                      ##  mas
e_dRA =  np.sqrt(e_RAg**2+e_RAi**2)                           ##  mas
e_dDE =  np.sqrt(e_DEg**2+e_DEi**2)                           ##  mas
C = (Cg*e_RAg*e_DEg + Ci*e_RAi*e_DEi)/(e_dRA*e_dDE)

#print 'Now write '
#fou = open('../data/qsocom.dat', 'w')
#print>>fou, \
#'''## Combined data of icrf2 catatlog and gaia qso.dat
### suffix 'g' for gaia data and 'i' for icrf2
### IERS  RA_g  DE_g  eRA*_g  eDE_g  RA_i  DE_i  eRA*_i  eDE_i  d_RA*  d_DE  e_dRA*  e_DE  D/N  FrameFix
###       "     "     mas     mas    "     "     mas     mas    mas    mas   mas     mas     ''' 
#
#for i in range(len(ICRFg)):
#    if ICRFg[i] in icrf2d:
#        flag = 'D'
#    else:
#        flag = 'N'
#        
#    print>>fou, ICRFg[i], '%14.6f  '*12%(RAg[i], DEg[i], e_RAg[i], e_DEg[i], \
#    RAi[i], DEi[i], e_RAi[i], e_DEi[i], d_RA[i], d_DE[i], e_dRA[i], e_dDE[i]),\
#    flag, '   %d'%Frot[i]
#    
#fou.close()

### some subsets
## All sources
inda = list(range(ICRFg.size))
## defining sources
indd = [ k for k in range(ICRFg.size) if ICRFg[k] in icrf2d ] 
## non-defining sources
indn = [ k for k in range(ICRFg.size) if ICRFg[k] not in icrf2d ] 
## Frot = 3
indr3= np.where(Frot==3)[0]
## Frot = 3
indr1= np.where(Frot==1)[0]

Ind = [inda, indd, indn, indr3]
Lab = ['All', 'defining', 'non-defining', 'fixed-frame']
## Notice!
# Units for RA and DE are arcsec and sould be translated to deg
#SpinFit(d_RA, d_DE, e_dRA, e_dDE,\
#                                np.deg2rad(RAg/3600),np.deg2rad(DEg/3600), flog)
#pmRA = d_RA
#pmDE = d_DE
#RA = np.deg2rad(RAg/3600)
#DE = np.deg2rad(DEg/3600)
RA = RAg/3600
DE = DEg/3600
#for i in range(len(Ind)):
#    print>>flog, 'subset %s: '%Lab[i]
#    ind0 = Ind[i]
#
#    d_RA0  = np.take(d_RA,  ind0)
#    d_DE0  = np.take(d_DE,  ind0)
#    e_dRA0 = np.take(e_dRA, ind0)
#    e_dDE0 = np.take(e_dDE, ind0)
#    RA0    = np.take(RA,    ind0)
#    DE0    = np.take(DE,    ind0)
#    C0     = np.take(C,     ind0)

## without considering the correlation between \alpha and \delta.    
#    print>>flog, 'Weighted: Normal.'
#    w, sig, corrcoef = RotationFit(d_RA0, d_DE0, e_dRA0, e_dDE0, RA0, DE0)

## considering the correlation between \alpha and \delta.
#    print>>flog, 'Weighted: Correlated.'
#    w, sig, corrcoef = RotationFitC(d_RA0, d_DE0, e_dRA0, e_dDE0, RA0, DE0, C0)
#    
#    [wx,   wy,  wz] = w
#    [ewx, ewy, ewz] = sig
#
#    print>>flog, 'the orientation/spin are(mas):\n :\n', \
#        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
#    print>>flog, '   correlation coefficients are:\n', corrcoef
    
### consider a possible bias in declination.    
#    print>>flog, 'Weighted: Normal. Consider a bias in declination.'
#    w, sig, corrcoef = RotationFitB(d_RA0, d_DE0, e_dRA0, e_dDE0, RA0, DE0)
#    [wx,   wy,  wz, b]  = w
#    [ewx, ewy, ewz, eb] = sig
#
#    print>>flog, 'the orientation are(mas):\n :\n', \
#        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
#    print>>flog, 'the bias in declination(mas):\n :\n', \
#        '    &$%.3f \pm %.3f$'%(b, eb)    
#    print>>flog, '   correlation coefficients are:\n', corrcoef

## fixed frame ones    
print('subset 260 fixed frame: ', file=flog)

d_RA0  = np.hstack((np.take( d_RA, indr3), np.take( d_RA, indr1)))
d_DE0  = np.hstack((np.take( d_DE, indr3), np.zeros(indr1.size)))
e_dRA0 = np.hstack((np.take(e_dRA, indr3), np.take(e_dRA, indr1)))
e_dDE0 = np.hstack((np.take(e_dDE, indr3), np.take(e_dDE, indr1)))
RA0    = np.hstack((np.take(   RA, indr3), np.take(   RA, indr1)))
DE0    = np.hstack((np.take(   DE, indr3), np.take(   DE, indr1)))
C0     = np.hstack((np.take(    C, indr3), np.take(    C, indr1)))

## without considering the correlation between \alpha and \delta.   
#print>>flog, 'Weighted: Normal.'
#w, sig, corrcoef = RotationFit(d_RA0, d_DE0, e_dRA0, e_dDE0, RA0, DE0)
## considering the correlation between \alpha and \delta.
#print>>flog, 'Weighted: Correlated.'
#w, sig, corrcoef = RotationFitC(d_RA0, d_DE0, e_dRA0, e_dDE0, RA0, DE0, C0)
#    
#[wx,   wy,  wz] = w
#[ewx, ewy, ewz] = sig
#
#print>>flog, 'the orientation/spin are(mas*yr-1):\n :\n', \
#        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
#print>>flog, '   correlation coefficients are:\n', corrcoef

### consider a possible bias in declination.    
#print>>flog, 'Weighted: Normal. Consider a bias in declination.'
#w, sig, corrcoef = RotationFitB(d_RA0, d_DE0, e_dRA0, e_dDE0, RA0, DE0)
#[wx,   wy,  wz, b]  = w
#[ewx, ewy, ewz, eb] = sig
#
#print>>flog, 'the orientation are(mas):\n :\n', \
#        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
#print>>flog, 'the bias in declination(mas):\n :\n', \
#        '    &$%.3f \pm %.3f$'%(b, eb)    
#print>>flog, '   correlation coefficients are:\n', corrcoef
##                
##### Mean and Std for offset.
##MdRA = np.average(d_RA, weights=e_dRA**-2)
##StdRA= np.std(d_RA/e_dRA**2)
##print>>flog, 'Mean and Std of d_RA are: %7.4fmas   %7.4fmas'%(MdRA, StdRA)
##
##MdDE = np.average(d_DE, weights=e_dDE**-2)
##StdDE= np.std(d_DE/e_dDE**2)
##print>>flog, 'Mean and Std of d_DE are: %7.4fmas   %7.4fmas'%(MdDE, StdDE)
##
#print>>flog, time.strftime('%Y-%m-%d %H:%M:%S Finished!',time.localtime(time.time()))
#flog.close()
#print 'Done!'
plt.plot(DE0, d_DE0, '.')
plt.ylim([-1.0, 1.0])
plt.show()
