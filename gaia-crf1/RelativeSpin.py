# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 23:44:44 2016

@author: Neo

Read the merged fits file of Tycho-2 and Tags catalog.

"""

from ReadMeg import readMeg
import time
from RotationFit import MeanWei

## Log file.
flog = open('../logs/spin.log', 'a')

### Tycho-2 stars
print>>flog, '''
################################################################################
### Spin between TGAS and Hip2 (All) begins.'''
### Spin between TGAS and Tycho-2 (All) begins.'''
print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))

#fname = '../data/merge_all.fits'
fname = '../data/merge_Hip2.fits'
tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
            plx, plx_err, l, b, g_mag = readMeg(fname, 'hip')
#            'tycho2-id')

dpmra = pmraG - pmraT
dpmde = pmdecG - pmdecT
ra = raG
de = decG

MeanWei(dpmra, dpmde, ra, de, flog)

### For M-K giants
print>>flog, '''
################################################################################
### Spin between TGAS and Hip2 (M-K giants only) begins.'''
### Spin between TGAS and Tycho-2 (M-K giants only) begins.'''
print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
#fname = '../data/merge_all_mk.fits'
fname = '../data/merge_Hip2_mk.fits'
tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
            plx, plx_err, l, b, g_mag = readMeg(fname, 'hip')
#            'tycho2-id')

dpmra = pmraG - pmraT
dpmde = pmdecG - pmdecT
ra = raG
de = decG

MeanWei(dpmra, dpmde, ra, de, flog)

### For O-B5 stars
print>>flog, '''
################################################################################
### Spin between TGAS and Hip2 (O-B5 stars only) begins.'''
### Spin between TGAS and Tycho-2 (O-B5 stars only) begins.'''
print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
#fname = '../data/merge_ob.fits'
fname = '../data/merge_Hip2_ob.fits'
tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
            plx, plx_err, l, b, g_mag = readMeg(fname, 'hip')
#            'tycho2-id')

dpmra = pmraG - pmraT
dpmde = pmdecG - pmdecT
ra = raG
de = decG

MeanWei(dpmra, dpmde, ra, de, flog)