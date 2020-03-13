# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 21:21:41 2016

@author: Neo

a new script to calculate orientation and spin

x, y, z -> wx, wy, wz/ex, ey, ez

"""

import numpy as np
sin = np.sin
cos = np.cos

def Cal(dalp, ddet, erra, errd, alp, det):
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
