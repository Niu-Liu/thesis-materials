# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:43:55 2015

Plot coodinate time series for radio sources. 
@author: Neo
"""
import numpy as np
import matplotlib.pyplot as plt
from fun import ADepoA, ADepoS
cos = np.cos

dat_dir = '../data/opa/'
res_dir = '../plot/timeseries/'

t0 = 2000.0

def tsplot(soun, pmra, pmdec, ra0, dec0):
    epo, ra, dec, era, edec = np.loadtxt(dat_dir+soun +'.dat', usecols=list(range(5)), unpack=True)
    if epo.size>1:
        epo = ADepoA(epo)
    else:
        epo = ADepoS(epo)
    if ra0 == 0.0:
        ra0 = ra[-1]
        dec0= dec[-1] 
    x, y1, err1, y2, err2 = epo, (ra-ra0)*3.6e6*cos(np.deg2rad(dec)), era, (dec-dec0)*3.6e6, edec
    x0 = t0 
    x1 = np.arange(1979.0, 2017.0, 0.1)
##  time series plot 
    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    ax0.errorbar(x, y1, yerr=err1, fmt='bo', markersize=3)    
    ax1.errorbar(x, y2, yerr=err2, fmt='bo', markersize=3)
###  for data points >=9:
    if pmra != 0.0:
        y3 = pmra*(x1-x0)/1.0e3
        y4 = pmdec*(x1-x0)/1.0e3
        ax0.plot(x1, y3, 'r')
        ax1.plot(x1, y4, 'r')
    
##  some details.    
    ax0.set_ylabel('R.A.(mas)')
    ax0.set_ylim([-50, 50])
    ax0.set_xlim([1979,2017])
    ax0.set_title(soun)
    ax1.set_ylabel('Dec(mas)')
    ax1.set_ylim([-50, 50])
#    plt.show()
    plt.savefig(res_dir+soun+'.eps', dpi=100)
    plt.close()
    
#tsplot('0434-188')
    
## read catalog file to get name of sources.
cat = '../list/opa.list'
soun = np.loadtxt(cat, dtype=str)
##  linear drift data.
apm = '../results/opa_all.apm'
pmRA, pmDE, RA0, DE0 = np.loadtxt(apm, usecols=(2,3,7,8), unpack=True)

for i in range(len(soun)):
    sou_name = soun[i]
    pmra, pmdec, ra0, dec0 = pmRA[i], pmDE[i], RA0[i], DE0[i]
## plot
    tsplot(sou_name, pmra, pmdec, ra0, dec0)
    
print('Done!')