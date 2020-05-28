# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:26:58 2016

@author: Neo
"""

from astropy.io import fits

def readMeg(fname, tagn='tycho2-id'):
        hdulist = fits.open(fname, memmap=True)
        tbdat = hdulist[1].data
#        print tbdat, hdulist[1].data
        
        tag = tbdat.field(tagn)
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
        

        return tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
            plx, plx_err, l, b, g_mag