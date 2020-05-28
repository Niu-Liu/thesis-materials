# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:04:22 2016

@author: Neo
"""

from astropy.io import fits
import numpy as np
import time

def extract_merge():
    '''
    extract O-B5 and M-K III data from Merged data.
    '''
    dat_dir = '/Users/Neo/Astron/Work/201609_HCRF/data/'
    
    flog = open('../logs/extract_merge.log', 'a')
    print>>flog, ''' extract O-B5 and M-K III data from Merged data.'''
    print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    
    tagmk = np.loadtxt(dat_dir+'tycho2_III231.mk', dtype=str, usecols=(0,))
    tagob = np.loadtxt(dat_dir+'tycho2_III231.ob', dtype=str, usecols=(0,))
    
    hdulist = fits.open(dat_dir+'merge_all.fits', memmap=True)
#    hdulist = fits.open(dat_dir+'merge_nd_all.fits', memmap=True)
    tbdat = hdulist[1].data
    tag0 = tbdat.field('tycho2-id')

    nummk = 0
    numob = 0
    N = len(tag0)    
 
    print 'Now it begins! '

## corresponding indice of MK giants and O-B5 stars.
    indmk = []
    indob = []

    for i in range(N):
        if tag0[i] in tagob:
            indob.append(i)
            numob += 1
        elif tag0[i] in tagmk:
            indmk.append(i)
            nummk += 1

        print '#: ' + '%d/%d '%(i,N)+ ' is completed.'
        
    print>>flog, '%d '%numob, ' O-B5 stars are found.'
    print>>flog, '%d '%nummk, ' MK giants are found.'

## read the merged data.
    raT  = tbdat.field('RA1')
    ra_errT = tbdat.field('e_RA1')
    decT = tbdat.field('DE1')
    dec_errT = tbdat.field('e_DE1')
    pmraT = tbdat.field('pmRA1')
    pmra_errT = tbdat.field('e_pmRA1')
    pmdecT = tbdat.field('pmDE1')
    pmdec_errT = tbdat.field('e_pmDE1')
        
    raG  = tbdat.field('RA2')
    ra_errG = tbdat.field('e_RA2')
    decG = tbdat.field('DE2')
    dec_errG = tbdat.field('e_DE2')
    pmraG = tbdat.field('pmRA2')
    pmra_errG = tbdat.field('e_pmRA2')
    pmdecG = tbdat.field('pmDE2')
    pmdec_errG = tbdat.field('e_pmDE2')
        
    plx = tbdat.field('Plx')
    plx_err = tbdat.field('e_Plx')  
    l = tbdat.field('GLON')
    b = tbdat.field('GLAT')
    g_mag = tbdat.field('phot_g_mean_mag')

## substract the common part from Tycho-2 
    print 'Now write data into fits format.'
    ind0 = [indob, indmk]
    typ0 = ['ob', 'mk']
    for i in range(len(typ0)):
        ind = ind0[i]
        typ = typ0[i]
#        print ind
        
        ra1 = np.take(raT, ind)
        dec1 = np.take(decT, ind)
        pmra1 = np.take(pmraT, ind)
        pmdec1 = np.take(pmdecT, ind)
        ra_err1 = np.take(ra_errT, ind)
        dec_err1 = np.take(dec_errT, ind)
        pmra_err1 = np.take(pmra_errT, ind)
        pmdec_err1 = np.take(pmdec_errT, ind)
        
        tagcom = np.take(tag0, ind)    
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
                fits.Column(name='Plx', format='F6.3', array=plx2),\
                fits.Column(name='e_Plx', format='F4.2', array=plx_err2),\
                fits.Column(name='GLON', format='F14.10', array=l2),\
                fits.Column(name='GLAT', format='F14.10', array=b2),\
                fits.Column(name='phot_g_mean_mag', format='D', array=g_mag2),\
                ])
    
## header information of FITS      
        prihdr = fits.Header()
        prihdr['Creator'] = 'Niu Liu'
        prihdr['COMMENT'] = "Parameters with suffix: 1 from Tycho-2 catalog, 2 from TGAS catalog"
        prihdr['COMMENT'] = "metadata please see /I/337 in cds"
        prihdu = fits.PrimaryHDU(header=prihdr)
        thdulist = fits.HDUList([prihdu, tbhdu])
#        thdulist.writeto(dat_dir+'merge_all_'+typ+'.fits') 
        thdulist.writeto(dat_dir+'merge_nd_all_'+typ+'.fits')

    print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    flog.close()

    print 'Done!'
    
extract_merge()