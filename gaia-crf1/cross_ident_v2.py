# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 22:18:51 2016

@author: Neo

corss identification between Revised Hipparcos and TGAS.
For every TGAS entries, we try to find its corresponding enteries in Hipparcos.
And classify them into three groups:
    All;
    K-M giants;
    O-B5 stars.
"""

import numpy as np
from astropy.io import fits
import time

dat_dir = '/Users/Neo/Astronomy/Work/201609_GDR1_CRF/data/'

def crossident_TT():
## log file.    
    flog = open('../logs/cross_ident.log', 'a')
    print>>flog, ''' Cross identification between TGAS and HIP2 begins.'''
    print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
## Hip2 data.
    fhip2 = dat_dir + 'hip2_single.dat'
    hip1  = np.loadtxt(fhip2, dtype=int, usecols=(0,))
## Gaia data.
    ftgas = dat_dir + 'TagsSource.hip'
    hip2  = np.loadtxt(ftgas, dtype=int, usecols=(0,))

## corresponding indice between catalog1 and catalog2
    indice1 = []
    indice2 = []

    num = 0
    N = len(hip2) 

    print 'cross idenfication begins! '
    for i in range(N):
        j = np.where(hip1==hip2[i])[0]
        if len(j):
            indice2.append(i)
            indice1.append(j[0])
            num += 1

        print '#: ' + '%d/%d '%(i,N)+ ' is completed.'
        
    print '%d '%num, 'common entries are found.'

    ra3, dec3, plx3, pmra3, pmdec3, ra_err3, dec_err3, plx_err3, pmra_err3, pmdec_err3 = \
        np.loadtxt(fhip2, unpack=True, usecols=range(1,11,1))
    ra4, ra_err4, dec4, dec_err4, plx4, plx_err4, pmra4, pmra_err4, pmdec4, pmdec_err4 = \
        np.loadtxt(ftgas, unpack=True, usecols=range(1,11,1))
         
    l, b, g_mag = np.loadtxt(ftgas, unpack=True, usecols=(21,22,23))

## substract the common part from Tycho-2    
    ra1 = np.take(ra3, indice1)
    dec1 = np.take(dec3, indice1)
    pmra1 = np.take(pmra3, indice1)
    pmdec1 = np.take(pmdec3, indice1)
    plx1 = np.take(plx3, indice1)
    ra_err1 = np.take(ra_err3, indice1)
    dec_err1 = np.take(dec_err3, indice1)
    plx_err1 = np.take(plx_err3, indice1)
    pmra_err1 = np.take(pmra_err3, indice1)
    pmdec_err1 = np.take(pmdec_err3, indice1)
    
    tagcom = np.take(hip2, indice2)    
    ra2 = np.take(ra4, indice2)
    dec2 = np.take(dec4, indice2)
    plx2 = np.take(plx4, indice2)
    pmra2 = np.take(pmra4, indice2)
    pmdec2 = np.take(pmdec4, indice2)
    ra_err2 = np.take(ra_err4, indice2)
    dec_err2 = np.take(dec_err4, indice2)
    plx_err2 = np.take(plx_err4, indice2)
    pmra_err2 = np.take(pmra_err4, indice2)
    pmdec_err2 = np.take(pmdec_err4, indice2)

    l = np.take(l, indice2)
    b = np.take(b, indice2)
    g_mag = np.take(g_mag, indice2)
    
    print 'Now write data into fits format.'

    tag = 'hip'
    tagfor = 'I6'
        
    tbhdu = fits.BinTableHDU.from_columns(\
            [fits.Column(name=tag, format=tagfor, array=tagcom), \
## data from Hip2 catalog.
            fits.Column(name='RA1', format='F14.10', array=ra1),\
            fits.Column(name='DE1', format='F14.10', array=dec1),\
            fits.Column(name='pmRA1', format='F9.3', array=pmra1),\
            fits.Column(name='pmDE1', format='F9.3', array=pmdec1),\
            fits.Column(name='e_RA1', format='F6.3', array=ra_err1),\
            fits.Column(name='e_DE1', format='F6.3', array=dec_err1),\
            fits.Column(name='e_pmRA1', format='F6.3', array=pmra_err1),\
            fits.Column(name='e_pmDE1', format='F6.3', array=pmdec_err1),\
#            fits.Column(name='Plx1', format='F6.3', array=plx1),\
#            fits.Column(name='e_Plx1', format='F4.2', array=plx_err1),\
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
            fits.Column(name='GLON', format='F14.10', array=l),\
            fits.Column(name='GLAT', format='F14.10', array=b),\
            fits.Column(name='phot_g_mean_mag', format='D', array=g_mag),\
            ])

## header information of FITS      
    prihdr = fits.Header()
    prihdr['Creator'] = 'Niu Liu'
    prihdr['COMMENT'] = "Parameters with suffix: 1 from Hip2 catalog, 2 from TGAS catalog"
    prihdr['COMMENT'] = "metadata please see /I/337 in cds"
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(dat_dir+'merge_Hip2.fits')        

    print>>flog, ': %d/%d entries are found. Completed!'%(num, N)
    print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    flog.close()

    print 'Done!'

crossident_TT()