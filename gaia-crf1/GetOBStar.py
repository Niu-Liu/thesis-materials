# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:44:36 2016

@author: Neo

get young O-B5 stars(age>10^7yrs)
"""

import numpy as np
from ReadMeg import readMeg
from astropy.io import fits

dat_dir = '/Users/Neo/Astron/Work/201609_HCRF/data/'
hiplist = ['hip_7.30.dat', 'hip_7.47.dat', 'hip_7.78.dat']
tyclist = ['tyc_7.30.dat', 'tyc_7.47.dat', 'tyc_7.78.dat']

## Mapping between of Hipparcos and Tycho-2 identifiers. 
TYC = np.loadtxt(dat_dir+'Hip_Tyc.map', dtype=str, usecols=(1,))
Hip = np.loadtxt(dat_dir+'Hip_Tyc.map', dtype=str, usecols=(0,))

## Array to store the Tycho identifiers
tyc0 = np.array([],dtype=str)

## read data in turns.
for i in range(len(hiplist)):
    temp = []
    hip = np.loadtxt(dat_dir+hiplist[i], dtype=str, usecols=(0,))
    
    for j in range(len(hip)):
        k = np.where(Hip==hip[j])[0]
        if len(k):
            temp.append(TYC[k[0]])
        else:
            print hip[j], " doesn't have a corresponding Tycho identifier! "
            
    tyc0 = np.hstack((tyc0, temp))
    
for i in range(len(tyclist)):
    tyc = np.loadtxt(dat_dir+tyclist[i], dtype=str, usecols=(0,))
    tyc0 = np.hstack((tyc0, tyc))
    
print len(tyc0)

fname = '../data/merge_all_ob.fits'
tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
            plx, plx_err, l, b, g_mag = readMeg(fname, 'tycho2-id')
            
ind =[]
for i in range(len(tag)):
    if tag[i] not in tyc0:
        ind.append(i)
            
ra1 = np.take(raT, ind)
dec1 = np.take(decT, ind)
pmra1 = np.take(pmraT, ind)
pmdec1 = np.take(pmdecT, ind)
ra_err1 = np.take(ra_errT, ind)
dec_err1 = np.take(dec_errT, ind)
pmra_err1 = np.take(pmra_errT, ind)
pmdec_err1 = np.take(pmdec_errT, ind)
        
tagcom = np.take(tag, ind)    
ra2 = np.take(raG, ind)
dec2 = np.take(decG, ind)
pmra2 = np.take(pmraG, ind)
pmdec2 = np.take(pmdecG, ind)
ra_err2 = np.take(ra_errG, ind)
dec_err2 = np.take(dec_errG, ind)
pmra_err2 = np.take(pmra_errG, ind)
pmdec_err2 = np.take(pmdec_errG, ind)
    
plx2 = np.take(plx, ind)
plx_err2 = np.take(plx_err, ind)
l2 = np.take(l, ind)
b2 = np.take(b, ind)
g_mag2 = np.take(g_mag, ind)

taglabel = 'tycho2-id'
tagfor = '12A'
            
tbhdu = fits.BinTableHDU.from_columns(\
                [fits.Column(name=taglabel, format=tagfor, array=tagcom), \
    ## data from Tycho-2 catalog.
                fits.Column(name='RA1', format='F14.10', array=ra1),\
                fits.Column(name='DE1', format='F14.10', array=dec1),\
                fits.Column(name='pmRA1', format='F9.3', array=pmra1),\
                fits.Column(name='pmDE1', format='F9.3', array=pmdec1),\
                fits.Column(name='e_RA1', format='F6.3', array=ra_err1),\
                fits.Column(name='e_DE1', format='F6.3', array=dec_err1),\
                fits.Column(name='e_pmRA1', format='F6.3', array=pmra_err1),\
                fits.Column(name='e_pmDE1', format='F6.3', array=pmdec_err1),\
## data from TGAS catalog.
                fits.Column(name='RA2', format='F14.10', array=ra2),\
                fits.Column(name='DE2', format='F14.10', array=dec2),\
                fits.Column(name='pmRA2', format='F9.3', array=pmra2),\
                fits.Column(name='pmDE2', format='F9.3', array=pmdec2),\
                fits.Column(name='e_RA2', format='F6.3', array=ra_err2),\
                fits.Column(name='e_DE2', format='F6.3', array=dec_err2),\
                fits.Column(name='e_pmRA2', format='F6.3', array=pmra_err2),\
                fits.Column(name='e_pmDE2', format='F6.3', array=pmdec_err2),\
                fits.Column(name='Plx', format='F6.3', array=plx2),\
                fits.Column(name='e_Plx', format='F4.2', array=plx_err2),\
                fits.Column(name='GLON', format='F14.10', array=l2),\
                fits.Column(name='GLAT', format='F14.10', array=b2),\
                fits.Column(name='phot_g_mean_mag', format='D', array=g_mag2),\
                ])
    
## header information of FITS      
prihdr = fits.Header()
prihdr['Creator'] = 'Niu Liu'
prihdr['COMMENT'] = " Data f O-B5 stars with age lager than 10^7yrs. "
prihdr['COMMENT'] = "Parameters with suffix: 1 from Tycho-2 catalog, 2 from TGAS catalog"
prihdr['COMMENT'] = "metadata please see /I/337 in cds"
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto(dat_dir+'merge_ob.fits')        

print 'Done!'