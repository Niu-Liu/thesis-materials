# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 20:39:13 2016

@author: Neo

structure index pick.

"""

import numpy as np
import time

dat_dir = '../data/'
##ICRF2 defining sources data file
icrf2f = 'icrf2-defining.dat'
defn = np.loadtxt(dat_dir+icrf2f, dtype=str, skiprows=20, usecols=(2,))
##input
inf = '../list/new_candidate1.cat'
soun = np.loadtxt(inf, dtype=str)
## output
cadf = '../list/new_candidate2.cat'
exf =  '../list/ex_def2.cat'
fcad = open(cadf,'w')
fexf = open(exf, 'w')
print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
    time.localtime(time.time())), file=fcad)
print('## Structure Index <=3', file=fcad)
print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
    time.localtime(time.time())), file=fexf)
print('## Structure Index =4', file=fexf)

## read structure index data
fil2 = '../data/BVID_SX_SI.dat'
l2 = np.loadtxt(fil2, dtype=str, usecols=(0,))
si = np.loadtxt(fil2, dtype=int, usecols=(1,))

## find the Sources with SI<=3
for i in range(soun.size):
    if soun[i] in l2:
        j = np.where(l2==soun[i])[0][0]
        if si[j]<=3:
                print(soun[i], file=fcad)
        elif soun[i] in defn:
                print(soun[i], file=fexf)
    else:
        print(soun[i], file=fcad)

fcad.close()
fexf.close()        