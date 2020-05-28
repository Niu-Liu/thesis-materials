# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 23:01:14 2016

@author: Neo

parameter fitting, using Oort-Lindblad equation:

k*u_l^* = (S1*sin(l) - S2*cos(l)/r + B*cos(b) + A*cos(2l)*cos(b)

"""

import numpy as np
sin = np.sin
cos = np.cos
k = 4.7407

def ParFit(pmlon, err, l, b, r):
    '''
    pmls = pml*cos(b), mas/yr
    r: kpc
    '''
    
    parS1 =  sin(l)/r
    parS2 = -cos(l)/r
    parA  =  cos(2*l)*cos(b)
    parB  =  cos(b)

    wgt = np.diag(err**-2)
    
    FT = np.vstack((parS1, parS2, parA, parB))
    F = np.transpose(FT)
    
## a*x = y
    a = np.dot(np.dot(FT, wgt), F)
    y = np.dot(np.dot(FT, wgt), k*np.transpose(pmlon))
    
    x = np.linalg.solve(a, y)    
    cov = np.linalg.inv(a)  
    sig = np.sqrt(cov.diagonal())
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(x)) for i in range(len(x))])
    corrcoef.resize((len(x), len(x)))
    
    return x, sig, corrcoef
    
#def Elim(pmRAs,pmDEs, e_pmRA,e_pmDE, n):
#    a = np.abs(pmRAs) - n*e_pmRA
#    b = np.abs(pmDEs) - n*e_pmDE
#    indice1 = np.where(a<0)
#    indice2 = np.where(b<0)
#    indice  = np.intersect1d(indice1, indice2)
##    pmRA = np.take(pmRA, indice)
##    pmDE = np.take(pmDE, indice)
##    e_pmRA = np.take(e_pmRA, indice)
##    e_pmDE = np.take(e_pmDE, indice)
##    RA = np.take(RA, indice)
##    DE = np.take(DE, indice)    
#
#    return indice
#    
#def loop(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE, flog):
#    w, sig, corrcoef = RotationFit(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE)
#
#    print>>flog, 'the orientation/spin are(mas/mas*yr-1):\n ',w
#    print>>flog, '   sigma are(mas/mas*yr-1):\n', sig
#    print>>flog, '   correlation coefficients are:\n', corrcoef
###  Calculate the residuals
#    wx = w[0]
#    wy = w[1]
#    wz = w[2]
#    
#    pmRAf = -sin(DE)*cos(RA)*wx-sin(DE)*sin(RA)*wy+cos(DE)*wz
#    pmDEf =  sin(RA)*wx-cos(RA)*wy
#    
### residuals
#    pmRAs = pmRAf - pmRA
#    pmDEs = pmDEf - pmDE    
#    
### generate the new array after eliminating outliers 
#    n = 2.6
#    indice = Elim(pmRAs,pmDEs, e_pmRA,e_pmDE, n)
#    npmRA = np.take(pmRA, indice)
#    npmDE = np.take(pmDE, indice)
#    ne_pmRA = np.take(e_pmRA, indice)
#    ne_pmDE = np.take(e_pmDE, indice)
#    nRA = np.take(RA, indice)
#    nDE = np.take(DE, indice) 
#    
#    print>>flog, 'Eliminating(%2.1f-sigma): %d Outliers.'%(n,len(pmRA)-len(npmRA))
#    
#    w, sig, corrcoef = RotationFit(npmRA,npmDE,ne_pmRA,ne_pmDE,nRA,nDE)

###  Calculate the residuals again
#    wx = w[0]
#    wy = w[1]
#    wz = w[2]
#    
#    npmRAf = -sin(nDE)*cos(nRA)*wx-sin(nDE)*sin(nRA)*wy+cos(nDE)*wz
#    npmDEf =  sin(nRA)*wx-cos(nRA)*wy
#    
### residuals
#    npmRAs = npmRAf - npmRA
#    npmDEs = npmDEf - npmDE
#
##    res = np.hstack((npmRAs/ne_pmRA**2, npmDEs/ne_pmDE**2))
#    wrmsRA = np.sqrt(np.sum(npmRAs**2/ne_pmRA**2)/np.sum(ne_pmRA**-2))
#    wrmsDE = np.sqrt(np.sum(npmDEs**2/ne_pmDE**2)/np.sum(ne_pmDE**-2))  
#    
#    print>>flog, 'After eliminating outliers, the orientation/spin are(mas/mas*yr-1):\n ',w
#    print>>flog, '   sigma are(mas/mas*yr-1):\n', sig
#    print>>flog, '   correlation coefficients are:\n', corrcoef
#    print>>flog, '   w.r.m.s(mas/mas*yr-1) is:\n', wrmsRA, wrmsDE
#    
#    return w, sig, corrcoef, wrmsRA, wrmsDE
#    
#def SpinFit(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE, flog):
####  First Loop, calculate the mean error without weights.
#    np.ones(len(pmRA))  
#    print>>flog, 'First, no weight apply to the data. Calculate the w.r.m.s'    
#    w, sig, corrcoef, wrmsRA, wrmsDE = loop(pmRA,pmDE,np.ones(len(pmRA)),np.ones(len(pmDE)),RA,DE, flog)
#    merrRA, merrDE = wrmsRA, wrmsDE
#    
###  Next, apply the weigths of M_err**2/err_i**2
#    print>>flog, 'Next, apply the weigths of M_err**2/err_i**2'
#    w, sig, corrcoef, wrmsRA, wrmsDE = loop(pmRA,pmDE,e_pmRA/merrRA,e_pmDE/merrDE,RA,DE, flog)

###############################################################################
## Main function        