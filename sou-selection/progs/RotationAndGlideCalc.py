# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 16:17:59 2016

Consider Rotation and Glide vectors simultaneously.

@author: Neo

Oct 25 2016: updated by Niu.
"""

import numpy as np
cos = np.cos
sin = np.sin
import time
from RotationAndGlideFitting import MeanWeiGliAndRot

## some path variables
res_dir = '../results/'
dat_fil = ['OV.rank', 'GR.rank']
res = ['OVrag.dat', 'GRrag.dat']

for d in range(len(dat_fil)):
    fout = open(res_dir + res[d],'w')
    print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
                               time.localtime(time.time())), file=fout)
    print('## Data file: %s'%dat_fil[d], file=fout)
    print('''## Global Spin and Glide vectors of selected source subsets.
##data format
##Num   w_tol   w_x   w_y   w_z
##            (Expectation | Std)
##unit: uas/yr
############################################################################''', file=fout)

    Sou = np.loadtxt(res_dir + dat_fil[d], dtype=str, usecols=(0,)) 
    V, V_err, RA, V_RA, Va_E, Dec, V_Dec, Vd_E = \
        np.loadtxt(res_dir + dat_fil[d], usecols=list(range(1,9)), unpack=True) 
## Unit: deg -> rad
    RA = np.deg2rad(RA)
    Dec = np.deg2rad(Dec)
    
    N = len(Sou)
    
    for i in range(100, N+1):
        sig = np.ones(i)
        ra, dec = RA[:i], Dec[:i]
        ua, ud = V_RA[:i], V_Dec[:i]
        
        x, sig, corrcoef = MeanWeiGliAndRot(ua, ud, ra, dec)
#        eua, ede = Va_E[:i], Vd_E[:i]
#        w, sig, corrcoef = RotationFit(ua, ud, eua, ede, ra, dec)
        w, g = x.reshape((2,3))
        sigw, sigg = sig.reshape((2,3))
        
        wx, wy, wz = w
        sigwx, sigwy, sigwz = sigw
        w = np.sqrt(sum(w**2))
        sigw = np.sqrt(sum(sigw**2))
        
        gx, gy, gz = g
        siggx, siggy, siggz = sigg
        g = np.sqrt(sum(g**2))
        sigg = np.sqrt(sum(sigg**2))
    
        print(("%3d"+"    %8.4f"*16)%\
            (i, w, sigw, wx, sigwx, wy, sigwy, wz, sigwz,\
                g, sigg, gx, siggx, gy, siggy, gz, siggz ), file=fout)
        
    fout.close()
        
print('Done!')    
