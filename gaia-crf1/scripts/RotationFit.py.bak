# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 22:56:22 2016

@author: Neo

Rotation Componets fitting.

Fitting Equation:
    d_pmRA = -w_x*sin(DE)*cos(RA) - w_y*sin(DE)*sin(RA) + w_z*cos(DE)
    d_pmDE = +w_x*sin(RA)         - w_y*cos(RA)

"""

import numpy as np
sin = np.sin
cos = np.cos
pi = np.pi

def RotCal(dalp, ddet, erra, errd, alp, det):
    '''
    unit: rad
    '''
    
    A11 = np.sum( sin(det)**2*cos(alp)**2/erra**2) + \
        np.sum( sin(alp)**2/errd**2)
    A12 = np.sum( sin(det)**2*cos(alp)*sin(alp)/erra**2) + \
        np.sum(-sin(alp)* cos(alp)/errd**2)
    A13 = np.sum(-sin(det)*cos(det)*cos(alp)/erra**2) + \
        0.0
    A22 = np.sum( sin(det)**2*sin(alp)**2/erra**2) + \
        np.sum( cos(alp)**2/errd**2)
    A23 = np.sum(-sin(det)*cos(det)*sin(alp)/erra**2) + \
        0.0
    A33 = np.sum(cos(det)**2/erra**2) + \
        0.0
        
    b1 = np.sum(-sin(det)*cos(alp)*dalp/erra**2) + \
        np.sum( sin(alp)*ddet/errd**2)
    b2 = np.sum(-sin(det)*sin(alp)*dalp/erra**2) + \
        np.sum(-cos(alp)*ddet/errd**2)
    b3 = np.sum( cos(det)*dalp/erra**2) + \
        0.0
        
    A = np.mat([[A11, A12, A13],
                [A12, A22, A23],
                [A13, A23, A33]])
                
    b = np.array([b1, b2, b3])
    
    w = np.linalg.solve(A, b)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    sig = np.array(sig).reshape(3,)
    
#    print type(sig)
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
                 
    return w, sig, corrcoef


def RotationFit(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE):
## Notice !!!
## unit for RA and DE are rad.
## partial of d_pmRA, d_pmDE with respect to w_x, w_y, w_z 
    N = len(pmRA)    
    parx1 = -sin(DE)*cos(RA)
    parx2 =  sin(RA)
    pary1 = -sin(DE)*sin(RA)
    pary2 = -cos(RA)
    parz1 =  cos(DE)
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

def RotationFitB(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE):
## A bias in decclination is considered.
    N = len(pmRA) 
    
    parx1 = -sin(DE)*cos(RA)
    parx2 =  sin(RA)
    pary1 = -sin(DE)*sin(RA)
    pary2 = -cos(RA)
    parz1 =  cos(DE)
    parz2 =  np.zeros(N)
    parb1 =  np.zeros(N)
    parb2 =  np.ones(N)
    
    wei1 = np.diag(e_pmRA**-2)
    wei2 = np.diag(e_pmDE**-2)
    
    F1T = np.vstack((parx1, pary1, parz1, parb1))
    F1  = np.transpose(F1T)
    
    F2T = np.vstack((parx2, pary2, parz2, parb2))
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
    
def RotationFitC(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE,C):
## Notice !!!
## This subroute takes the correlation between \alpha and \delta into consideration.
    N = len(pmRA)    
    parx1 = -sin(DE)*cos(RA)
    parx2 =  sin(RA)
    pary1 = -sin(DE)*sin(RA)
    pary2 = -cos(RA)
    parz1 =  cos(DE)
    parz2 =  np.zeros(N)
    
    wei1 = np.diag(e_pmRA**-2)
    wei2 = np.diag(e_pmDE**-2)
    
    weic = np.diag(C/e_pmRA/e_pmDE)
    
    F1T = np.vstack((parx1, pary1, parz1))
    F1  = np.transpose(F1T)
    
    F2T = np.vstack((parx2, pary2, parz2))
    F2  = np.transpose(F2T)

    a1 = np.dot(np.dot(F1T, wei1), F1)
    b1 = np.dot(np.dot(F1T, wei1), pmRA)
    
    a2 = np.dot(np.dot(F2T, wei2), F2)
    b2 = np.dot(np.dot(F2T, wei2), pmDE)

## coupling terms
    c1 = np.dot(np.dot(F1T, weic), F2)
    d1 = np.dot(np.dot(F1T, weic), pmDE)
    
    c2 = np.dot(np.dot(F2T, weic), F1)
    d2 = np.dot(np.dot(F2T, weic), pmRA)
    
    A = a1 + a2 + c1 + c2
    b = b1 + b2 + d1 + d2

    w = np.linalg.solve(A, b)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
                 
    return w, sig, corrcoef
    
def CalcRes(pmRA, pmDE, RA, DE, w, sig):
    '''
    Calculate the residuals.
    '''
    [wx,   wy,  wz] = w
    [ewx, ewy, ewz] = sig
    
##  Calculate the residuals
    pmRAf = -sin(DE)*cos(RA)*wx-sin(DE)*sin(RA)*wy+cos(DE)*wz
    pmDEf =  sin(RA)*wx-cos(RA)*wy
        
## residuals
    pmRAs = pmRAf - pmRA
    pmDEs = pmDEf - pmDE 
    varRA = np.sum(pmRAs**2)/(len(pmRAs)-1)
    varDE = np.sum(pmDEs**2)/(len(pmDEs)-1)
    wrms = np.sqrt((varRA + varDE)/2)
    
    return pmRAs, pmDEs, wrms
    
def Elim(pmRAs,pmDEs, e_pmRA,e_pmDE, n):  
    a = np.abs(pmRAs) - n*e_pmRA
    b = np.abs(pmDEs) - n*e_pmDE
    indice1 = np.where(a<0)
    indice2 = np.where(b<0)
    indice  = np.intersect1d(indice1, indice2) 

    return indice

def MeanWei(pmRA,pmDE, RA,DE, flog):
    Num1 = 1   # number of outliers in the last iteration
    Num2 = 0   # number of outliers in this iteration
    cout = 0   # number of iterations
    n = 2.6    # n-sigma filter
    wrms = 1   # mean uncertainty
    
    npmRA, npmDE, nRA, nDE = pmRA, pmDE, RA, DE
        
    while Num1 != Num2:
        Num1 = Num2
        w, sig, corrcoef = \
            RotCal(npmRA,npmDE,wrms*np.ones(len(npmRA)),wrms*np.ones(len(npmDE)),nRA,nDE)
        pmRAs, pmDEs, wrms = CalcRes(pmRA, pmDE, RA, DE, w, sig) 
## generate the new array after eliminating outliers 
        indice = Elim(pmRAs,pmDEs, wrms*np.ones(len(pmRAs)),wrms*np.ones(len(pmRAs)), n)
        npmRA  = np.take(pmRA, indice)
        npmDE  = np.take(pmDE, indice)
        nRA    = np.take(RA, indice)
        nDE    = np.take(DE, indice) 
        Num2   = len(pmRA)-len(npmRA)
        cout  += 1
        print>>flog, 'Iteration %d: %d'%(cout, Num2)
        
    print>>flog, 'Eliminating(%2.1f-sigma): %d/%d (objects/Outliers).'%(n,len(pmRA), Num2)  
    
    [wx,   wy,  wz] = w
    [ewx, ewy, ewz] = sig

    print>>flog, '''
    After elimination: '''
    print>>flog, 'the orientation/spin are(mas*yr-1):\n :\n', \
        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
    print>>flog, '   correlation coefficients are:\n', corrcoef
    
def NormalWei(pmRA, pmDE, RA, DE, pmRAerr, pmDEerr, flog):
    Num1 = 1   # number of outliers in the last iteration
    Num2 = 0   # number of outliers in this iteration
    cout = 0   # number of iterations
    n = 2.6    # n-sigma filter
#    wrms = 1   # mean uncertainty
    
    npmRA, npmDE, nRA, nDE, npmRAerr, npmDEerr = pmRA, pmDE, RA, DE, pmRAerr, pmDEerr
        
    while Num1 != Num2:
        Num1 = Num2
        w, sig, corrcoef = \
            RotCal(npmRA,npmDE,npmRAerr,npmDEerr,nRA,nDE)
        pmRAs, pmDEs, wrms = CalcRes(pmRA, pmDE, RA, DE, w, sig) 
## generate the new array after eliminating outliers 
        indice = Elim(pmRAs,pmDEs, pmRAerr, pmDEerr, n)
        npmRA  = np.take(pmRA, indice)
        npmDE  = np.take(pmDE, indice)
        nRA    = np.take(RA, indice)
        nDE    = np.take(DE, indice) 
        npmRAerr  = np.take(pmRAerr, indice)
        npmDEerr  = np.take(pmDEerr, indice)
        Num2   = len(pmRA)-len(npmRA)
        cout  += 1
        print>>flog, 'Iteration %d: %d'%(cout, Num2)
        
    print>>flog, 'Eliminating(%2.1f-sigma): %d/%d (objects/Outliers).'%(n,len(pmRA), Num2)  
    
    [wx,   wy,  wz] = w
    [ewx, ewy, ewz] = sig

    print>>flog, '''
    After elimination: '''
    print>>flog, 'the orientation/spin are(mas*yr-1):\n :\n', \
        '    &$%.3f \pm %.3f$'*3%(wx, ewx, wy, ewy, wz, ewz)
    print>>flog, '   correlation coefficients are:\n', corrcoef