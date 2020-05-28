# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:39:38 2016

According the rank list based on linear drift, calculate the rotation of three axis

@author: Neo

Oct 21 2016: updated by Niu.
"""

import numpy as np
cos = np.cos
sin = np.sin
import time
from RotationFit import MeanWeiRotation, RotationFit

## some path variables
res_dir = '../results/'
dat_fil = ['OV.rank', 'GR.rank']
#dat_fil = ['OV1.rank', 'GR1.rank']
res = ['OVrot.dat', 'GRrot.dat']

for d in range(len(dat_fil)):
    fout = open(res_dir + res[d],'w')
    print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
                               time.localtime(time.time())), file=fout)
    print('## Data file: %s'%dat_fil[d], file=fout)
    print('''## Global Spin Vector of selected source subsets.
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
    
    N= len(Sou)
    
    for i in range(100,len(Sou)+1):
        sig = np.ones(i)
        ra, dec = RA[:i], Dec[:i]
        ua, ud = V_RA[:i], V_Dec[:i]
        
        w, sig, corrcoef = MeanWeiRotation(ua, ud, ra, dec)
#        eua, ede = Va_E[:i], Vd_E[:i]
#        w, sig, corrcoef = RotationFit(ua, ud, eua, ede, ra, dec)
        wx, wy, wz = w
        sigwx, sigwy, sigwz = sig
        w = np.sqrt(sum(w**2))
        sigw = np.sqrt(sum(sig**2))
    
        print(("%3d"+"    %8.4f"*8)%(i, w, sigw, wx, sigwx, wy, sigwy, wz, sigwz), file=fout)
        
fout.close()
print('Done!')