# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:22:54 2016

pre-selection(criterion 3)

@author: neo

Oct 17 2016, Niu: update.
"""

import numpy as np
import time

##ICRF2 defining sources data file
icrf2f = '../data/icrf2-defining.dat'
defn = np.loadtxt(icrf2f, dtype=str, skiprows=20, usecols=(2,))

in_fil = '../results/variance.dat'
soun = np.loadtxt(in_fil, dtype=str, usecols=(0,))
Stda, Asta, WASa, Stdd, AStd, WASd = \
    np.loadtxt(in_fil, usecols=list(range(3,6))+list(range(7,10)), unpack=True) 

out_fil = '../list/new_candidate3.cat'
fou = open(out_fil,'w')
ex_fil = '../list/ex_def3.cat'
fex = open(ex_fil, 'w')

print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
    time.localtime(time.time())), file=fou)
print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
    time.localtime(time.time())), file=fex)
print('## ICRF2 defining sources, ejected.', file=fex)

for i in range(len(soun)):
    if max(Stda[i],Stdd[i],WASa[i],WASd[i]) < 10:
        print(soun[i], file=fou)
    elif soun[i] in defn:
        print(soun[i], '  %8.2f'*4%(Stda[i],Stdd[i],WASa[i],WASd[i]), file=fex)
        
fou.close()
fex.close()

print('Done!')