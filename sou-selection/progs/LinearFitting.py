#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:17:17 2016

@author: Neo

linear least square fitting, a specially weighting approach used.

"""

import numpy as np

t0 = 2000.0
Mpm = 1.0e3
MEpm= 1.0e3

def mean(x, sig):
    N = x.size
## Standard deviation of residual series.
    if N == 1:
        ave, std_ave = x[0], sig[0]
    else:    
        p = sig**-2
## weighted average
        ave = sum(p*x)/sum(p)
## residual series.
        v = np.sqrt(p)*(x-ave)
        std_ave = np.sqrt(sum(v**2)/sum(p)/(N-1))
    
    return [ave, std_ave]

def linlsq(x, y, sig):
## Fitting curves: y = a*x + b 
    para = x
    parb = np.ones(x.size)
    
    Ft = np.vstack((para, parb))
    P = np.diag(sig**-2)
    F = np.transpose(Ft)
    
## normal equation:  A*p = B
    A = np.dot(np.dot(Ft, P), F)
    B = np.dot(np.dot(Ft, P), y)
    
    p = np.linalg.solve(A, B)   
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    sig = np.array(sig).reshape(len(p),)
    
    a, b = p
    siga, sigb = sig
    
    return [a, b, siga, sigb]

def Elim(RA, DE, e_RA,e_DE, n):
    a = np.abs(RA) - n*e_RA
    b = np.abs(DE) - n*e_DE
    indice1 = np.where(a<0)
    indice2 = np.where(b<0)
    indice  = np.intersect1d(indice1, indice2)
    return indice
    
def pmFit(epo, ra, era, dec, edec):
    t = epo - t0

    if t.size:
        if t.size < 3:
            pmra, pmdec = 0.0, 0.0
#            pmra, pmdec = Mpm, Mpm
            epmra, epmdec = MEpm, MEpm
            ra0,  era0 = mean(ra,  era)
            dec0, edec0= mean(dec,edec)
        else:
##  first fitting.
            [pmra,  ra0,  epmra,  era0 ] = linlsq(t, ra,  era)
            [pmdec, dec0, epmdec, edec0] = linlsq(t, dec, edec)
##  calculate the residuals
            vra = ra - (pmra*t + ra0)
            vdec= dec- (pmdec*t+dec0)
##  eliminate outliers: 3.0 sigma
            ind = Elim(vra, vdec, era, edec, 3.0)
##  substract the good point from original observational data.
            nt  = np.take(t,  ind)
            nra = np.take(ra, ind)
            ndec= np.take(dec,ind)
            nera = np.take(era, ind)
            nedec= np.take(edec,ind)
            
            if nt.size:
                if nt.size < 3:
                    pmra, pmdec = 0.0, 0.0
#                    pmra, pmdec = Mpm, Mpm
                    epmra, epmdec = MEpm, MEpm
                    ra0,  era0 = mean(nra,  nera)
                    dec0, edec0= mean(ndec,nedec)
                else:
                    [pmra,  ra0,  epmra,  era0 ] = linlsq(nt, nra,  nera)
                    [pmdec, dec0, epmdec, edec0] = linlsq(nt, ndec,nedec)
            else:
                pmra, pmdec = 0.0, 0.0
#                pmra, pmdec = Mpm, Mpm
                epmra, epmdec = MEpm, MEpm
                ra0, dec0 = 0.0, 0.0
                era0, edec0 = MEpm, MEpm
    else:
        pmra, pmdec = 0.0, 0.0
#        pmra, pmdec = Mpm, Mpm
        epmra, epmdec = MEpm, MEpm
        ra0, dec0 = 0.0, 0.0
        era0, edec0 = MEpm, MEpm
        
    return [pmra,  ra0,  epmra,  era0, pmdec, dec0, epmdec, edec0]
    
def APM(epo0, ra0, era0, dec0, edec0):
##  Only data points over 1979.0-2015.0 are used
    ind = np.where(epo0<=2015.0)[0]
    epo = np.take(epo0,  ind)
    ra  = np.take(ra0, ind)
    dec = np.take(dec0,ind)
    era = np.take(era0, ind)
    edec= np.take(edec0,ind)
##  we devide the series into 3 samples:
##  1979.0~1990.0, 1990.0~2009.0, 2009.0~2015.0
    if epo.size <= 20:
#        pmRA, pmDE = 0.0, 0.0
        if epo.size >= 3:
            [pmRA, RA0, epmRA, eRA0] = linlsq(epo, ra,  era)
            [pmDE, DE0, epmDE, eDE0] = linlsq(epo, dec,edec)
        else:
            pmRA, pmDE = 0.0, 0.0
            epmRA, epmDE = MEpm, MEpm
            RA0, eRA0 = mean(ra0,  era0 )
            DE0, eDE0 = mean(dec0, edec0)
    else:
##  And if the data point is too far away from its neighborhood, it will be ejected
#        ind = [ i for i in range(epo.size-1) \
#            if max([np.fabs(ra[i]-ra[i-1]), np.fabs(ra[i]-ra[i+1]), \
#                    np.fabs(dec[i]-dec[i-1]), np.fabs(dec[i]-dec[i+1])])<1000 ]
#        epo  = np.take(epo,  ind)
#        ra = np.take(ra, ind)
#        dec= np.take(dec,ind)
#        era = np.take(era, ind)
#        edec= np.take(edec,ind)
        
        t1, t2 = 1990.0, 2009.0
        ind1 = np.where(epo < t1)[0]
    
        ind20= np.where(epo >=t1)[0]
        ind21= np.where(epo < t2)[0]
        ind2 = np.intersect1d(ind20, ind21)
        
        ind3 = np.where(epo>=t2)[0]
    
        ind = [ind1, ind2, ind3]
## fitting the apparent proper motion seperatively.
        pmRA, pmDE, RA0, DE0, epmRA, epmDE, eRA0, eDE0  = np.zeros(8)

        
        for i in range(len(ind)):
            epos = np.take(epo, ind[i])
            ras  = np.take(ra,  ind[i])
            decs = np.take(dec, ind[i])
            eras = np.take(era, ind[i])
            edecs= np.take(edec,ind[i])
            
            epmra, epmdec, era0, edec0  = np.ones(4)*MEpm
            [pmra,  ra0,  epmra,  era0, pmdec, dec0, epmdec, edec0] = \
                    pmFit(epos, ras, eras, decs, edecs)
            
#            print pmra,  ra0,  epmra,  era0, pmdec, dec0, epmdec, edec0
            pmRA += pmra/epmra**2
            pmDE += pmdec/epmdec**2
            
            epmRA+= 1.0/epmra**2
            epmDE+= 1.0/epmdec**2
                
            RA0 += ra0/era0**2
            DE0 += dec0/edec0**2
            eRA0+= 1.0/era0**2
            eDE0+= 1.0/edec0**2
            
        if RA0 and DE0 and epmRA and epmDE:    
#            print pmRA, epmRA            
            pmRA = pmRA/epmRA
            pmDE = pmDE/epmDE
            epmRA = 1.0/np.sqrt(epmRA)
            epmDE = 1.0/np.sqrt(epmDE)
                    
            RA0 = RA0/eRA0
            DE0 = DE0/eDE0
            eRA0 = 1.0/np.sqrt(eRA0)
            eDE0 = 1.0/np.sqrt(eDE0)
            
        else:
            pmRA, pmDE = 0.0, 0.0
#            pmRA, pmDE = 5.0, 5.0
            epmRA, epmDE = MEpm, MEpm
            RA0, eRA0 = mean(ra,  era )
            DE0, eDE0 = mean(dec, edec)
            
#        earray = np.array([epmRA, epmDE, eRA0, eDE0])
#        epmRA, epmDE, eRA0, eDE0 = np.where(earray>1.0e6, 1.0e6, earray)
    
    return [pmRA, pmDE, RA0, DE0, epmRA, epmDE, eRA0, eDE0]