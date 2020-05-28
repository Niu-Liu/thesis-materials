# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 22:18:51 2016

@author: Neo

corss identification between Tycho-2 and TGAS.
For every TGAS entries, we try to find its corresponding enteries in Tycho-2.

Sep 18 2016, LN: modified crsid_TT, trying to make it effctive
"""

import numpy as np
from astropy.io import fits
import time

dat_dir = '/Users/Neo/Astron/Work/201609_HCRF/data/'

def crsid_TT(s):
    fty2 = dat_dir + 'tyc2.'+s
    ftgas= dat_dir + 'TagsSource.'+s
    
    flog = open('../logs/cross_ident.log', 'a')
    print>>flog, ''' Cross identification between TGAS and Tycho-2 begins.'''
    print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    
    tag1 = np.loadtxt(fty2, dtype=str, usecols=(0,))
    tag2 = np.loadtxt(ftgas, dtype=str, usecols=(0,))
#    flag = list(tag1)
    
#    ra3, dec3, pmra3, pmdec3, ra_err3, dec_err3, pmra_err3, pmdec_err3 = \
#    np.loadtxt(fty2, usecols=(range(1,9)), unpack=True)
#    
#    ra4, ra_err4, dec4, dec_err4, plx, plx_err, pmra4, pmra_err4, pmdec4, pmdec_err4, l, b = np.loadtxt(ftgas, usecols=(range(1,11)+[21,22]), unpack=True)
#                ra_dec, ra_para, ra_pmra, ra_pmdec, \
#                dec_prl, dec_pmra, dec_pmdec, \
#                prl_pmra, prl_pmdec, pmra_pmdec, \
#    
#    N = len(tag2)
#    ra1 = np.zeros(N)
#    dec1 = np.zeros(N)
#    pmra1 = np.zeros(N)
#    pmdec1 = np.zeros(N)
#    ra_err1 = np.zeros(N)
#    dec_err1 = np.zeros(N)
#    pmra_err1 = np.zeros(N)
#    pmdec_err1 = np.zeros(N)
#    
#
#    ra2 = np.zeros(N)
#    ra_err2 = np.zeros(N)
#    dec2 = np.zeros(N)
#    dec_err2 = np.zeros(N)
#    plx = np.zeros(N)
#    plx_err = np.zeros(N)
#    pmra2 = np.zeros(N)
#    pmra_err2 = np.zeros(N)
#    pmdec2 = np.zeros(N)
#    pmdec_err2 = np.zeros(N)
#    l = np.zeros(N)
#    b  = np.zeros(N)

    num = 0
    N = len(tag1)    
 
    print 'cross idenfication begins! '

#    for i in range(N):
#        if tag2[i] in flag:
#            j = flag.index(tag2[i])
#            ra1[num], dec1[num], pmra1[num], pmdec1[num], \
#                    ra_err1[num], dec_err1[num], pmra_err1[num], pmdec_err1[num] =\
#            ra3[j], dec3[j], pmra3[j], pmdec3[j], \
#                    ra_err3[j], dec_err3[j], pmra_err3[j], pmdec_err3[j]
#            
#            ra2[num], ra_err2[num], dec2[num], dec_err2[num], \
#                    plx[num], plx_err[num], \
#                    pmra2[num], pmra_err2[num], pmdec2[num], pmdec_err2[num], l[num], b[num] = \
#            ra4[i], ra_err4[i], dec4[i], dec_err4[i], \
#            plx[i], plx_err[i], \
#            pmra4[i], pmra_err4[i], pmdec4[i], pmdec_err4[i], l[i], b[i]
#            num += 1
#
#        print s + ': ' + '%f '%(i*1.000/N)+ '% is completed.'


## corresponding indice between catalog1 and catalog2
    indice1 = []
    indice2 = []

    for i in range(N):
#    for i in range():
        j = np.where(tag2==tag1[i])[0]
        if len(j):
            indice1.append(i)
            indice2.append(j[0])
            num += 1

        print s + ': ' + '%d/%d '%(i,N)+ ' is completed.'
        
    print '%d '%num, 'common entries are found.'
    
#    print indice1[0], indice2[0]
        
    ra3, dec3, pmra3, pmdec3, ra_err3, dec_err3, pmra_err3, pmdec_err3 = \
    np.loadtxt(fty2, usecols=(range(1,9)), unpack=True)
    
    ra4, ra_err4, dec4, dec_err4, plx, plx_err, \
    pmra4, pmra_err4, pmdec4, pmdec_err4, l, b , g_mag\
    = np.loadtxt(ftgas, usecols=(range(1,11)+[21,22,23]), unpack=True)
#    print ra3[i]
#    print ra4[j]

## substract the common part from Tycho-2    
    ra1 = np.take(ra3, indice1)
    dec1 = np.take(dec3, indice1)
    pmra1 = np.take(pmra3, indice1)
    pmdec1 = np.take(pmdec3, indice1)
    ra_err1 = np.take(ra_err3, indice1)
    dec_err1 = np.take(dec_err3, indice1)
    pmra_err1 = np.take(pmra_err3, indice1)
    pmdec_err1 = np.take(pmdec_err3, indice1)
#    print ra1
    
    tagcom = np.take(tag2, indice2)    
    ra2 = np.take(ra4, indice2)
    dec2 = np.take(dec4, indice2)
    pmra2 = np.take(pmra4, indice2)
    pmdec2 = np.take(pmdec4, indice2)
    ra_err2 = np.take(ra_err4, indice2)
    dec_err2 = np.take(dec_err4, indice2)
    pmra_err2 = np.take(pmra_err4, indice2)
    pmdec_err2 = np.take(pmdec_err4, indice2)
#    print ra2

    plx = np.take(plx, indice2)
    plx_err = np.take(plx_err, indice2)
    l = np.take(l, indice2)
    b = np.take(b, indice2)
    g_mag = np.take(g_mag, indice2)
    
    print 'Now write data into fits format.'

## write the data into fits
    if s =='hip':
        tag = 'hip'
        tagfor = '6I'
    else:
        tag = 'tycho2-id'
        tagfor = '12A'
        
    tbhdu = fits.BinTableHDU.from_columns(\
            [fits.Column(name=tag, format=tagfor, array=tagcom), \
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
            fits.Column(name='Plx', format='F6.3', array=plx),\
            fits.Column(name='e_Plx', format='F4.2', array=plx_err),\
            fits.Column(name='GLON', format='F14.10', array=l),\
            fits.Column(name='GLAT', format='F14.10', array=b),\
            fits.Column(name='phot_g_mean_mag', format='D', array=g_mag),\
            ])

## header information of FITS      
    prihdr = fits.Header()
    prihdr['Creator'] = 'Niu Liu'
    prihdr['COMMENT'] = "Parameters with suffix: 1 from Tycho-2 catalog, 2 from TGAS catalog"
    prihdr['COMMENT'] = "metadata please see /I/337 in cds"
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(dat_dir+'merge_'+s+'.fits')        

    print>>flog, s+': %d/%d entries are found. Completed!'%(num, N)
    print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    flog.close()

    print 'Done!'

crsid_TT('hip')
crsid_TT('tyc')
