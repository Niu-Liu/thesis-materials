# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:49:28 2016

@author: Neo
"""

import numpy as np
import time

res_dir = '../results/'
dat_fil = ['OV.rank', 'GR.rank']
rot_fil = ['OVrot.dat', 'GRrot.dat']
sim_fil = ['OVsim.dat', 'GRsim.dat']
res = 'Selected.info'
#ICRF2 defining sources data file
icrf2f = '../data/icrf2-defining.dat'
defn = np.loadtxt(icrf2f, dtype=str, skiprows=20, usecols=(2,))
## data
Souo = np.loadtxt(res_dir + dat_fil[0], usecols=(0,), dtype=str)
Deco = np.loadtxt(res_dir + dat_fil[0], usecols=(6,))
Soug = np.loadtxt(res_dir + dat_fil[1], usecols=(0,), dtype=str)
Decg = np.loadtxt(res_dir + dat_fil[1], usecols=(6,))
 
w1, ew1, wx1, ewx1, wy1, ewy1, wz1, ewz1 = \
    np.loadtxt(res_dir + rot_fil[0], usecols=list(range(1,9)), unpack=True)
g1, eg1 = np.loadtxt(res_dir + sim_fil[0], usecols=(1,2), unpack=True)

w2, ew2, wx2, ewx2, wy2, ewy2, wz2, ewz2 = \
    np.loadtxt(res_dir + rot_fil[1], usecols=list(range(1,9)), unpack=True)
g2, eg2 = np.loadtxt(res_dir + sim_fil[1], usecols=(1,2), unpack=True)       
## For Our selected subsets.
#ind = np.array([223, 368, 194, 347])
ind = np.array([120, 223, 130, 194])
w  = [  w1[ind[0]],  w1[ind[1]],  w2[ind[2]],  w2[ind[3]] ]
ew = [ ew1[ind[0]], ew1[ind[1]], ew2[ind[2]], ew2[ind[3]] ]
wx = [ wx1[ind[0]], wx1[ind[1]], wx2[ind[2]], wx2[ind[3]] ]
ewx= [ewx1[ind[0]],ewx1[ind[1]],ewx2[ind[2]],ewx2[ind[3]] ]
wy = [ wy1[ind[0]], wy1[ind[1]], wy2[ind[2]], wy2[ind[3]] ]
ewy= [ewy1[ind[0]],ewy1[ind[1]],ewy2[ind[2]],ewy2[ind[3]] ]
wz = [ wz1[ind[0]], wz1[ind[1]], wz2[ind[2]], wz2[ind[3]] ]
ewz= [ewz1[ind[0]],ewz1[ind[1]],ewz2[ind[2]],ewz2[ind[3]] ]

g  = [  g1[ind[0]],  g1[ind[1]],  g2[ind[2]],  g2[ind[3]] ]
eg = [ eg1[ind[0]], eg1[ind[1]], eg2[ind[2]], eg2[ind[3]] ]

ind = ind + 100
Dec1 = Deco[:ind[0]+1] 
Dec2 = Deco[:ind[1]+1] 
Dec3 = Decg[:ind[2]+1] 
Dec4 = Decg[:ind[3]+1]
MeanDec = [np.mean(Dec1), np.mean(Dec2), np.mean(Dec3), np.mean(Dec4)]
Sou = [ Souo[:ind[0]+1], Souo[:ind[1]+1], Soug[:ind[2]+1], Soug[:ind[3]+1] ]           
com = [ [ j for j in range(ind[i]) if Sou[i][j] in defn] for i in range(ind.size) ]
N = [ len(com[i]) for i in range(ind.size) ]

## output
fout = open(res_dir + res,'w')
print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
            time.localtime(time.time())), file=fout)
print('''## information of 4 selected subsets.
##data format
##Subsets   w_tol   w_x   w_y   w_z  g_tol   MeanDec  Num_ICRF2
##            (Expectation | Std)
##unit: uas/yr
############################################################################''', file=fout)
for i in range(ind.size):
    print('%dSou'%ind[i] + \
        '   %4.1f \pm %4.1f'*5%(w[i], ew[i], wx[i], ewx[i], \
        wy[i], ewy[i], wz[i], ewz[i], g[i], eg[i]) + \
        '   %6.2f   %3d'%(MeanDec[i], N[i]), file=fout)

fout.close()
print('Done!')