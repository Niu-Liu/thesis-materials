# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 21:00:13 2016

@author: Neo
"""

import numpy as np
import matplotlib.pyplot as plt

res_dir = '../results/'
rot_fil = ['OVrot.dat', 'GRrot.dat']
sim_fil = ['OVsim.dat', 'GRsim.dat']
## For 4 subsets in previous works.
w3 = np.loadtxt(res_dir + 'SpecialSets.rot', usecols=(1,))
g3 = np.loadtxt(res_dir + 'SpecialSets.sim', usecols=(1,))

Num = np.loadtxt(res_dir + rot_fil[0], usecols=(0,), dtype=int)
##For the Overall Rank
w1 = np.loadtxt(res_dir + rot_fil[0], usecols=(1,))
g1 = np.loadtxt(res_dir + sim_fil[0], usecols=(1,))

w2 = np.loadtxt(res_dir + rot_fil[1], usecols=(1,))
g2 = np.loadtxt(res_dir + sim_fil[1], usecols=(1,))

wm = max(max(w1), max(w2), max(w3))
gm = max(max(g1), max(g2), max(g3))

## fraction of w and g in combinated q 
n1 = 0.5
n2 = 0.5
## Normalized and Calculate the combined q
q3 = w3/wm*n1 + g3/gm*n2
q1 = w1/wm*n1 + g1/gm*n2
q2 = w2/wm*n1 + g2/gm*n2

i = Num
y = np.ones(Num.size)

plt.plot(i, q1, 'b', linewidth = 3.0 )
plt.plot(i, q2, 'r', linewidth = 3.0 )
plt.plot(i, q3[0]*y, '--', label = '212 ICRF', linewidth = 5.0 )
plt.plot(i, q3[1]*y, '--', label = '295 ICRF', linewidth = 5.0 )
plt.plot(i, q3[2]*y, '--', label = '247 MFV', linewidth = 5.0 )
plt.plot(i, q3[3]*y, '--', label = '260 AMS', linewidth = 2.0 )
plt.xlim([100, max(i)])
plt.ylim([0, 1.0])
plt.ylabel('$Q$', fontsize = 20)
plt.xlabel('No. Sources', fontsize = 20)
ind = [120, 223, 130, 194]
#ind = [223, 368, 194, 347]
plt.plot(i[ind[0]], q1[ind[0]]-0.02, '^b', markersize=10)
plt.plot(i[ind[1]], q1[ind[1]]-0.02, '^b', markersize=10)
plt.plot(i[ind[2]], q2[ind[2]]-0.02, '^r', markersize=10)
plt.plot(i[ind[3]], q2[ind[3]]-0.02, '^r', markersize=10)
plt.legend()
plt.show()
plt.savefig('../plot/Quality.eps', dpi = 100)
plt.close()