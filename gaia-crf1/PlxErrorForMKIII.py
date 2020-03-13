#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 09:11:53 2016

@author: Neo

plx uncertainty statistics for M-K giants.
"""

import numpy as np
import time
from astropy.io import fits

## Log file.
flog = open('../logs/PlxStatistics.log', 'a')
print>>flog, time.strftime('## The scripts runs at %Y-%m-%d %H:%M:%S',\
                           time.localtime(time.time()))

## Load proper motion data of M-K giants.
fname = '../data/merge_all_mk.fits'
hdulist = fits.open(fname, memmap=True)
tbdat = hdulist[1].data  
plx = tbdat.field('Plx')
plx_err = tbdat.field('e_Plx') 
print>>flog, 'All Entries: %d'%len(plx)       
   
## First eject the object with a minus plx 
indp = np.where(plx>0)[0]
plxp = np.take(plx, indp)
plx_errp = np.take(plx_err, indp)
print>>flog, 'Entries with a minus plx: %d'%(len(plx)-len(plxp))

## plx_err/plx ratio statistics.
rat = plx_errp/plxp*100     ## percentages

print>>flog, '## err_plx/plx staistics:'
## short statistics
flt = np.arange(10, 60, 10)  ## filter for err_plx/plx
for i in range(flt.size):
    indf = np.where(rat<= flt[i])[0]
    print>>flog, '%2d : %5d(%3.1f)'%(flt[i], indf.size, indf.size*100.0/rat.size)

## A cumulative distribution of plx_err/plx ratio
import matplotlib.pyplot as plt
### Histogram
x = rat
num_bins = 5000000
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, facecolor='green', alpha=0.5,\
                            normed=1, histtype='step', cumulative=True)
plt.xlabel('$\sigma_{plx}/plx$(%)', fontsize=16)
#plt.ylabel('%', fontsize=18)
plt.yticks([0.2, 0.4, 0.6, 0.8, 0.862],\
               ['20', '40', '60', '80', '86.2'], fontsize=14)
plt.xlim([0, 100])
plt.xticks(range(0, 101, 10), fontsize=12)
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],\
               ['0', '', '20', '30', '40', '',  '60', '', '80', '', '100'], fontsize=14)
plt.vlines(30, 0, 0.862, linestyles='dashed', colors='r')
plt.hlines(0.862, 0, 30, linestyles='dashed', colors='r')
plt.plot(30, 0.862, 'ro')
plt.text(-4, 0.95, '%', fontsize=14)
### Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
#plt.show()
plt.savefig('../plot/PlxErr2Plx.eps', dpi=100)
plt.close()

### plx_err/plx ratio vs. plx
#plt.plot(rat, plxp, '.', markersize=1)
#plt.xlabel('$\sigma_{plx}/plx$(%)', fontsize=18)
#plt.ylabel('Plx(mas)', fontsize=18)
#plt.xlim([0,500])
#plt.ylim([0, 30])
#plt.show()
#plt.savefig('../plot/ErrRito2Plx.eps', dpi=100)
#plt.close()
print 'Done!'