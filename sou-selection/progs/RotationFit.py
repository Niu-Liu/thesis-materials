# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 22:56:22 2016

@author: Neo

Rotation Componets fitting.

Fitting Equation:
    d_pmRA = -w_x*sin(DE)*cos(RA) - w_y*sin(DE)*sin(RA) + w_z*cos(DE)
    d_pmDE = +w_x*sin(RA)         - w_y*cos(RA)
    
Oct 21 2016: updated by Niu

"""

import numpy as np
sin = np.sin
cos = np.cos
pi = np.pi

def RotationFit(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE):
## Notice !!!
## unit for RA and DE are rad.
## partial of d_pmRA, d_pmDE with respect to w_x, w_y, w_z 
    N = pmRA.size    
    parx1 = -sin(DE)*cos(RA)
    pary1 = -sin(DE)*sin(RA)
    parz1 =  cos(DE)
    parx2 =  sin(RA)
    pary2 = -cos(RA)    
    parz2 =  np.zeros(N)
    
    wei1 = np.diag(e_pmRA**-2)
    wei2 = np.diag(e_pmDE**-2)
    
    F1T = np.vstack((parx1, pary1, parz1))
    F1  = np.transpose(F1T)
    
    F2T = np.vstack((parx2, pary2, parz2))
    F2  = np.transpose(F2T)

    a1 = np.dot(np.dot(F1T, wei1), F1)
    b1 = np.dot(np.dot(F1T, wei1), pmRA)
    
    a2 = np.dot(np.dot(F2T, wei2), F2)
    b2= np.dot(np.dot(F2T, wei2), pmDE)
    
    A = a1 + a2
    b = b1 + b2

    w = np.linalg.solve(A, b)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
                 
    return w, sig, corrcoef
    
def Elim(pmRAs,pmDEs, e_pmRA,e_pmDE, n):  
    a = np.abs(pmRAs) - n*e_pmRA
    b = np.abs(pmDEs) - n*e_pmDE
    indice1 = np.where(a<0)
    indice2 = np.where(b<0)
    indice  = np.intersect1d(indice1, indice2)  

    return indice
    
#def NormalRotation(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE):
#    w, sig, corrcoef = RotationFit(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE)
###  Calculate the fitting values.
#    wx = w[0]
#    wy = w[1]
#    wz = w[2]       
#    pmRAf = -sin(DE)*cos(RA)*wx-sin(DE)*sin(RA)*wy+cos(DE)*wz
#    pmDEf =  sin(RA)*wx-cos(RA)*wy    
### residuals
#    pmRAs = pmRAf - pmRA
#    pmDEs = pmDEf - pmDE  
#    ## generate the new array after eliminating outliers 
#    n = 2.6
#    indice = Elim(pmRAs,pmDEs, e_pmRA,e_pmDE, n)
#    npmRA = np.take(pmRA, indice)
#    npmDE = np.take(pmDE, indice)
#    nRA = np.take(RA, indice)
#    nDE = np.take(DE, indice) 
#    ne_pmRA = np.take(e_pmRA, indice) 
#    ne_pmDE = np.take(e_pmDE, indice) 
#    w, sig, corrcoef = RotationFit(npmRA,npmDE,ne_pmRA,ne_pmDE,nRA,nDE)
#    return w, sig, corrcoef
   
def MeanWeiRotation(pmRA,pmDE, RA,DE):
###  First Loop, calculate the mean error without weights. 
    N =  pmRA.size 
    w, sig, corrcoef = RotationFit(pmRA,pmDE,np.ones(N),np.ones(N),RA,DE)
    ##  Calculate the residuals
    wx = w[0]
    wy = w[1]
    wz = w[2]
        
    pmRAf = -sin(DE)*cos(RA)*wx -sin(DE)*sin(RA)*wy +cos(DE)*wz
    pmDEf =  sin(RA)*wx         -cos(RA)*wy
        
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
    
    w, sig, corrcoef = RotationFit(npmRA,npmDE,np.ones(len(npmRA)),np.ones(len(npmDE)),nRA,nDE)
    ##  Calculate the residuals
    wx = w[0]
    wy = w[1]
    wz = w[2]
        
    pmRAf = -sin(DE)*cos(RA)*wx -sin(DE)*sin(RA)*wy +cos(DE)*wz
    pmDEf =  sin(RA)*wx         -cos(RA)*wy
        
    ## residuals
    pmRAs = pmRAf - pmRA
    pmDEs = pmDEf - pmDE  
    wrmsRA = np.sqrt(np.sum(pmRAs**2)/(len(pmRAs)-1))
    wrmsDE = np.sqrt(np.sum(pmDEs**2)/(len(pmRAs)-1))
    wrms = np.sqrt((wrmsRA**2 + wrmsDE**2)/2)
    merr = wrms
    
#    ##  Next, apply the weigths of 1/M_err**2
    w, sig, corrcoef = RotationFit(pmRA,pmDE,np.ones(len(pmRA))*merr,\
                                            np.ones(len(pmRA))*merr, RA, DE)
#    ##  Calculate the residuals
#    wx = w[0]
#    wy = w[1]
#    wz = w[2]
#        
#    pmRAf = -sin(DE)*cos(RA)*wx-sin(DE)*sin(RA)*wy+cos(DE)*wz
#    pmDEf =  sin(RA)*wx-cos(RA)*wy
#        
#    ## residuals
#    pmRAs = pmRAf - pmRA
#    pmDEs = pmDEf - pmDE           
#      
#    ## eliminating outliers again
#    n = 2.6
#    indice = Elim(pmRAs,pmDEs, np.ones(len(pmRA))*merr, np.ones(len(pmRA))*merr, n)
#    npmRA = np.take(pmRA, indice)
#    npmDE = np.take(pmDE, indice)
#    nRA = np.take(RA, indice)
#    nDE = np.take(DE, indice) 
        
#    w, sig, corrcoef = RotationFit(npmRA,npmDE,np.ones(len(npmRA))*merr,np.ones(len(npmDE))*merr,nRA,nDE)

    return w, sig, corrcoef