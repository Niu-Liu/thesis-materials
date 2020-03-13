#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:55:05 2016

@author: Neo

proper motion fitting.

"""

import numpy as np
sin = np.sin
cos = np.cos
k =4.74047

def PMfit(pml, pmb, pmerrl, pmerrb, l, b, r):
    '''
    Actually here 'pml' = pml*cos(b), which is usually labelled as \mu_{l*}
    
    A*x = b
    x = (X Y Z A B)^T, A: 5x5 matrix
    '''
    kpml = k*pml
    kpmb = k*pmb
    
    N = len(pml)  
    
    parl1 =  sin(l)/r
    parl2 = -cos(l)/r
    parl3 =  np.zeros(N)
    parl4 =  cos(2*l)*cos(b)
    parl5 =  cos(b)
    
    parb1 =  cos(l)*sin(b)/r
    parb2 =  sin(l)*sin(b)/r
    parb3 = -cos(b)/r
    parb4 = -0.5*sin(2*l)*sin(2*b)
    parb5 =  np.zeros(N)    
    
    wei1 = np.diag(pmerrl**-2)
    wei2 = np.diag(pmerrb**-2)
    
    M1T = np.vstack(\
        (parl1, parl2, parl3, parl4, parl5))
    M1  = np.transpose(M1T)
    
    M2T = np.vstack(\
        (parb1, parb2, parb3, parb4, parb5))
    M2  = np.transpose(M2T)

    a1 = np.dot(np.dot(M1T, wei1), M1)
    b1 = np.dot(np.dot(M1T, wei1), kpml)
    
    a2 = np.dot(np.dot(M2T, wei2), M2)
    b2 = np.dot(np.dot(M2T, wei2), kpmb)
    
    A = a1 + a2
    B = b1 + b2

    x = np.linalg.solve(A, B)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    sig = np.array(sig).reshape(len(x),)
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(x)) for i in range(len(x))])
    corrcoef.resize((len(x), len(x)))
                 
    return x, sig, corrcoef

def PMfit1(pml, pmb, pmerrl, pmerrb, l, b, r):
    '''
    Here we consider 9 parameters x = (S1, S2, S3, D_32, D_13, D_21, D+12, D+13, D+32)^T
    D_21, D+12 similar to Oort constant B, A.
    '''
    kpml = k*pml
    kpmb = k*pmb
    
    N = len(pml)  
    
    parl1 =  sin(l)/r
    parl2 = -cos(l)/r
    parl3 =  np.zeros(N)
    parl4 = -cos(l)*sin(b)
    parl5 = -sin(l)*sin(b)
    parl6 =  cos(b)
    parl7 =  cos(2*l)*cos(b)
    parl8 = -sin(l)*sin(b)
    parl9 =  cos(l)*sin(b)
    
    parb1 =  cos(l)*sin(b)/r
    parb2 =  sin(l)*sin(b)/r
    parb3 = -cos(b)/r
    parb4 =  sin(l)
    parb5 = -cos(l)
    parb6 =  np.zeros(N)
    parb7 = -0.5*sin(2*l)*sin(2*b)
    parb8 =  cos(l)*cos(2*b)
    parb9 =  sin(l)*cos(2*b)    
    
    wei1 = np.diag(pmerrl**-2)
    wei2 = np.diag(pmerrb**-2)
    
    M1T = np.vstack(\
        (parl1, parl2, parl3, parl4, parl5, parl6, parl7, parl8, parl9))
    M1  = np.transpose(M1T)
    
    M2T = np.vstack(\
        (parb1, parb2, parb3, parb4, parb5, parb6, parb7, parb8, parb9))
    M2  = np.transpose(M2T)

    a1 = np.dot(np.dot(M1T, wei1), M1)
    b1 = np.dot(np.dot(M1T, wei1), kpml)
    
    a2 = np.dot(np.dot(M2T, wei2), M2)
    b2 = np.dot(np.dot(M2T, wei2), kpmb)
    
    A = a1 + a2
    B = b1 + b2

    x = np.linalg.solve(A, B)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    sig = np.array(sig).reshape(len(x),)
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(x)) for i in range(len(x))])
    corrcoef.resize((len(x), len(x)))
                 
    return x, sig, corrcoef
