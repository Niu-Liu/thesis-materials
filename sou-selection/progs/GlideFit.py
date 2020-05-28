#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:50:47 2016

@author: Neo

Glide Fitting.

Fitting Equation:
    d_pmRA = +g_x*sin(RA)         - g_y*cos(RA)
    d_pmDE = +g_x*sin(DE)*cos(RA) + g_y*sin(DE)*sin(RA) - g_z*cos(DE)

"""

import numpy as np
sin = np.sin
cos = np.cos
pi  = np.pi

def GlideFit(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE):
## Notice !!!
## unit for RA and DE are rad.
## partial of d_pmRA, d_pmDE with respect to g_x, g_y, g_z 
    N = pmRA.size   
    
    parx3 =  sin(RA)
    pary3 = -cos(RA)
    parz3 =  np.zeros(N)
    
    parx4 = +sin(DE)*cos(RA)
    pary4 = +sin(DE)*sin(RA)
    parz4 = -cos(DE)
    
    
    wei1 = np.diag(e_pmRA**-2)
    wei2 = np.diag(e_pmDE**-2)
    
    F1T = np.vstack((parx3, pary3, parz3))
    F1  = np.transpose(F1T)
    
    F2T = np.vstack((parx4, pary4, parz4))
    F2  = np.transpose(F2T)

    a1 = np.dot(np.dot(F1T, wei1), F1)
    b1 = np.dot(np.dot(F1T, wei1), pmRA)
    
    a2 = np.dot(np.dot(F2T, wei2), F2)
    b2= np.dot(np.dot(F2T, wei2), pmDE)
    
    A = a1 + a2
    b = b1 + b2

    g = np.linalg.solve(A, b)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(g)) for i in range(len(g))])
    corrcoef.resize((len(g), len(g)))
                 
    return g, sig, corrcoef
    
def Elim(pmRAs,pmDEs, e_pmRA,e_pmDE, n):
    a = np.abs(pmRAs) - n*e_pmRA
    b = np.abs(pmDEs) - n*e_pmDE
    indice1 = np.where(a<0)
    indice2 = np.where(b<0)
    indice  = np.intersect1d(indice1, indice2)
    
    return indice
   
def MeanWeiGlide(pmRA,pmDE, RA,DE):
###  First Loop, calculate the mean error without weights.  
    g, sig, corrcoef = GlideFit(pmRA,pmDE,np.ones(len(pmRA)),np.ones(len(pmDE)),RA,DE)
##  Calculate the residuals
    gx = g[0]
    gy = g[1]
    gz = g[2]

    pmDEf =  sin(RA)*gx-cos(RA)*gy        
    pmRAf = -sin(DE)*cos(RA)*gx-sin(DE)*sin(RA)*gy+cos(DE)*gz
      
## residuals
    pmRAs = pmRAf - pmRA
    pmDEs = pmDEf - pmDE  
    wrmsRA = np.sqrt(np.sum(pmRAs**2)/(len(pmRAs)-1))
    wrmsDE = np.sqrt(np.sum(pmDEs**2)/(len(pmRAs)-1))
    wrms = np.sqrt((wrmsRA**2 + wrmsDE**2)/2)
    
## generate the new array after eliminating outliers 
    n = 2.6
    indice = Elim(pmRAs,pmDEs, wrms*np.ones(len(pmRAs)),wrms*np.ones(len(pmRAs)), n)
    npmRA = np.take(pmRA, indice)
    npmDE = np.take(pmDE, indice)
    nRA = np.take(RA, indice)
    nDE = np.take(DE, indice) 
    
    g, sig, corrcoef = GlideFit(npmRA,npmDE,np.ones(len(npmRA)),np.ones(len(npmDE)),nRA,nDE)
##  Calculate the residuals
    gx = g[0]
    gy = g[1]
    gz = g[2]
        
    pmDEf =  sin(RA)*gx-cos(RA)*gy        
    pmRAf = -sin(DE)*cos(RA)*gx-sin(DE)*sin(RA)*gy+cos(DE)*gz
        
    ## residuals
    pmRAs = pmRAf - pmRA
    pmDEs = pmDEf - pmDE  
    wrmsRA = np.sqrt(np.sum(pmRAs**2)/(len(pmRAs)-1))
    wrmsDE = np.sqrt(np.sum(pmDEs**2)/(len(pmRAs)-1))
    wrms = np.sqrt((wrmsRA**2 + wrmsDE**2)/2)
    merr = wrms    
        
##  Next, apply the weigths of 1/M_err**2
    g, sig, corrcoef = GlideFit(pmRA,pmDE,np.ones(len(pmRA))*merr,\
                                            np.ones(len(pmRA))*merr, RA, DE)
     
    return g, sig, corrcoef