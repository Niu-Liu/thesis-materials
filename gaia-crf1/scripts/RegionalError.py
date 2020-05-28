# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 22:45:50 2016

@author: Neo

regional differences

"""

import numpy as np
import matplotlib.pyplot as plt

## Nodes of RA and DE
sra = 10
sde = 10
Nra = np.arange([  0, 370, sra])    ## deg
Nde = np.arange([-90, 100, sde])    ## deg
lra = len(Nra)
lde = len(Nde)
## coordinates of center of grids
Mra = Nra[:lra-1] + sra/2
Mde = Nde[:lde-1] + sde/2

def compart(l1, l2):
    '''
    find the common part of l1 and l2, return a list of common part.
    '''
    com = []
    for i in range(len(l1)):
        if l1[i] in l2:
            com.append(l1[i])
    return com
    
def elim(xo, xc, err):
    '''
    eliminate outliers that statisfied the condition:
    abs(xo(obs) - xc(cal))> n*err
    xo, xc should be arrays.
    '''
    n =2.6
    res = np.abs(xo - xc) - n*err
    ind = np.where(res<0)
    
    return ind
    
def WgtMean(x, err):
    wgt = err*-2
    mean = np.sum(x*wgt)/np.sum(wgt)
    
    return mean 
    
def CalMean(x, err):
    if len(x)>3:
        mean = WgtMean(x, err)
        ind = elim(x, mean, err)
        xn = np.take(x, ind)
        errn = np.take(err, ind)
        mean = WgtMean(xn, errn)
    else:
        mean = 0.0
    return mean

def RegErr(RA, DE, x, y, e_x, e_y):
    '''
    compute the regional differences for a two-dimensional vector (x, y)
    '''
    Regx = np.zeros([lra-1, lde-1])
    Regy = np.zeros([lra-1, lde-1])
    
    for i in range(lra-1):
        id1 = np.where(RA>Nra[i])
        id2 = np.where(RA<Nra[i+1])
        idra = compart(id1, id2)
        
        for j in range(lde-1):
            id3 = np.where(DE>Nde[i])
            id4 = np.where(DE<Nde[i+1])
            idde = compart(id3, id4)
            
            idcom = compart(idra, idde)
            
            xp = np.take(x, idcom)
            yp = np.take(y, idcom)
            e_x = np.take(e_x, idcom)
            e_y = np.take(e_y, idcom)
            
#            mxp = CalMean(xp, e_x)
#            myp = CalMean(yp, e_y)
            
            Regx[i,j] = CalMean(xp, e_x)
            Regy[i,j] = CalMean(yp, e_y)
            
    return Regx, Regy
    
def RegPlot(Regx, Regy):
    plt.figure(figsize=(12,6))
    for i in range(lra-1):
            for j in range(lde-1):
                plt.arrow(Mra[i], Mde[i], Regx[i], Regy[j])
    
    plt.xlim([  0, 360])
    plt.xticks(Nra)
    plt.ylim([-90,  90])
    plt.xticks(Nde)
    
    plt.grid()
    plt.savefig('../plot/RegDiff.eps', dpi=100)