# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 09:12:26 2016

@author: Neo

Rotation and Glide fitting, first order

Fitting Equation:
    d_pmRA = -w_x*sin(DE)*cos(RA) - w_y*sin(DE)*sin(RA) + w_z*cos(DE) \
             +g_x*sin(RA)         - g_y*cos(RA)
    d_pmDE = +w_x*sin(RA)         - w_y*cos(RA)                       \
             +g_x*sin(DE)*cos(RA) + g_y*sin(DE)*sin(RA) - g_z*cos(DE)
    
"""

import numpy as np
sin = np.sin
cos = np.cos
pi = np.pi

def RGFit(pmRA,pmDE,e_pmRA,e_pmDE,RA,DE):
## Notice !!!
## unit for RA and DE are rad.
## partial of d_pmRA, d_pmDE with respect to w_x, w_y, w_z, g_x, g_y, g_z 
    N = len(pmRA)    
    parx1 = -sin(DE)*cos(RA)
    parx2 =  sin(RA)
    pary1 = -sin(DE)*sin(RA)
    pary2 = -cos(RA)
    parz1 =  cos(DE)
    parz2 =  np.zeros(N)
    
    parx4 = +sin(DE)*cos(RA)
    parx3 =  sin(RA)
    pary4 = +sin(DE)*sin(RA)
    pary3 = -cos(RA)
    parz4 = -cos(DE)
    parz3 =  np.zeros(N)
    
    wei1 = np.diag(e_pmRA**-2)
    wei2 = np.diag(e_pmDE**-2)
    
    F1T = np.vstack((parx1, pary1, parz1, parx3, pary3, parz3))
    F1  = np.transpose(F1T)
    
    F2T = np.vstack((parx2, pary2, parz2, parx4, pary4, parz4))
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