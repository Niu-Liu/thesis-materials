# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 15:55:03 2016

@author: Neo

Oct 7: ICRF2 defining sources in blue points while others in red ones.
"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.ticker import NullLocator
from matplotlib.lines import Line2D

pi = np.pi
cos = np.cos
sin = np.sin
sqrt = np.sqrt

## data
res_dir = '../results/'
dat_fil = ['OV.rank', 'GR.rank']
ind = np.array([223, 194]) + 100
RAo, DEo = np.loadtxt(res_dir + dat_fil[0], usecols=(3,6), unpack=True)
RAg, DEg = np.loadtxt(res_dir + dat_fil[1], usecols=(3,6), unpack=True)
SNo = np.loadtxt(res_dir + dat_fil[0], usecols=(0,), dtype=str)
SNg = np.loadtxt(res_dir + dat_fil[1], usecols=(0,), dtype=str)

RA0 = [RAo[:ind[0]+1], RAg[:ind[1]+1]] 
DE0 = [DEo[:ind[0]+1], DEg[:ind[1]+1]]
SN0 = [SNo[:ind[0]+1], SNg[:ind[1]+1]]

#ICRF2 defining sources data file
icrf2f = '../data/icrf2-defining.dat'
defn = np.loadtxt(icrf2f, dtype=str, skiprows=20, usecols=(2,))

#Hammer Projection
def projection(alpha, delta):
    x = []
    y = []
    
    for i in range(len(alpha)):
        ra = np.deg2rad(alpha[i])
        de = np.deg2rad(delta[i])
        
        B  = sqrt(1 + cos(de)*cos(ra/2))
        xi = 2*np.sqrt(2)*cos(de)*sin(ra/2)/B
        yi = np.sqrt(2)*sin(de)/B
        
        x.append(xi)
        y.append(yi)
        
    return [x,y]

def SphericPlot(RA, DE, Sou, figname):
    N = 8
    ra_s = 360.0/N
    de_s = 180.0/N
    
    ra0 = np.arange(-180, 180+ra_s, ra_s)
    de0 = np.arange( -90,  90+de_s, de_s)

    ra1 = np.arange(-180, 185, 5)
    de1 = np.arange( -90,  95, 5)

#plt.axes(figsize=(12, 8), frameon=False)
    plt.axes([0, 0, 1, 1], frameon=False)

##plot the grid
    for i in range(len(ra0)):
        ra2 = np.ones(len(de1))*ra0[i]
        de2 = de1
        [x0, y0] = projection(ra2, de2)
        plt.plot(x0, y0, 'k-', linewidth = 0.2)
        
    for i in range(len(de0)):
        ra2 = ra1
        de2 = np.ones(len(ra1))*de0[i]
        [x0, y0] = projection(ra2, de2)
        plt.plot(x0, y0, 'k-', linewidth = 0.2)        
    plt.xlim(-3.1, 3.1)
    plt.ylim(-1.6, 1.6)
    plt.text(-3.05, -0.03, ' 0h')
    plt.text( 2.85, -0.03, '24h')
    plt.text(-0.10, -1.50, '-90$^o$')
    plt.text(-0.10,  1.45, ' 90$^o$')
    plt.xticks([])
    plt.yticks([])
    
    [x, y] = projection(np.array(RA)-180, np.array(Dec))
    for i in range(Sou.size):
        if Sou[i] in defn:
            plt.plot(x[i], y[i], 'bo', fillstyle='none')
        else:
            plt.plot(x[i], y[i], 'bo')
    plt.savefig('../plot/'+figname+'SouDis.eps', dpi=100)
    plt.close()

for i in range(ind.size):
    RA, Dec, SN = RA0[i], DE0[i], SN0[i]
    SphericPlot(RA, Dec, SN, str(ind[i]))