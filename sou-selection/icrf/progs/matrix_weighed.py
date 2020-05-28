# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:52:47 2015

@author: Neo
"""

import numpy as np
import fun

def coefficence_matrix(x,y,sig):
        L = len(x)
        F = np.empty([L,2], dtype=float)
        Y = np.empty([L,1], dtype=float)
        P = np.zeros([L,L], dtype=float)
        
        for i in range(L):
            F[i,0] = x[i]
            F[i,1] = 1
            P[i,i] = sig[i]**-2
            Y[i,0] = y[i]
    
        N = (F.T.dot(P)).dot(F)
        b = (F.T.dot(P)).dot(Y)
        
        return [N,b]
        
def solve_matrix_equation(N,b):
#Solve the normal equation, return the eastimated parameters and errors.
    N = np.mat(N)
    b = np.mat(b)
    
    D = N.I
    x = D.dot(b)
    sig = [np.sqrt(D[i,i]) for i in range(len(D))]
#    print x
    return [x,sig]
        
def matrix_weighed_fitting(alp,siga,det,sigd,epo): 
    t0 = fun.ADepo(51545.0) 
# the origin of epoch, corresponding to JD2000.0
    L = len(epo)
    tNode = [1990.0, 2009.0]
# set nodes, dividing the whole data into groups,
# each considered as an independent catalog.
    ind = [0,0,0,L] 
#index of time series, corresponding to nodes.
    for i in range(len(tNode)):
        t0 = tNode[i]
        j = 0
        while (j < L and epo[j] <= t0):
            j += 1
        
        ind[i+1] = j
        
#    print ind
        
    N1 = np.zeros([2,2], dtype=float)
    b1 = np.zeros([2,1], dtype=float)
    
    N2 = np.zeros([2,2], dtype=float)
    b2 = np.zeros([2,1], dtype=float)
    
    for i in range(len(ind)-1):
        x = epo[ind[i]:ind[i+1]]
        if len(x) >= 3:
            y1 = alp[ind[i]:ind[i+1]]
            sig1 = siga[ind[i]:ind[i+1]]
            [N10,b10] = coefficence_matrix(x-t0,y1,sig1)
            
            y2 = det[ind[i]:ind[i+1]]
            sig2 = sigd[ind[i]:ind[i+1]]
            [N20,b20] = coefficence_matrix(x-t0,y2,sig2)
            
        else:
            N10 = np.zeros([2,2], dtype=float)
            b10 = np.zeros([2,1], dtype=float)
            
            N20 = np.zeros([2,2], dtype=float)
            b20 = np.zeros([2,1], dtype=float)
            
#        print N10,b10,N20,b20
        
        N1 += N10
        b1 += b10
        
        N2 += N20
        b2 += b20
        
#        print N1,b1,N2,b2
        
    [ra,sigra]=solve_matrix_equation(N1,b1)
    [de,sigde]=solve_matrix_equation(N2,b2)    
    
    return [ra,sigra,de,sigde]
