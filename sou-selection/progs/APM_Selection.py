# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 21:58:57 2016

@author: neo

Oct 19, 2016: updated by Niu.
"""

import numpy as np
import time

## generate the list of candidates and ejected ICRF2 defining sources.
out_fil = '../list/new_candidate4.cat'
fou = open(out_fil,'w')
ex_fil = '../list/ex_def4.cat'
fex = open(ex_fil, 'w')

print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
    time.localtime(time.time())), file=fou)
print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
    time.localtime(time.time())), file=fex)
print('## ICRF2 defining sources, ejected.', file=fex)

# APM data file.
apm = '../results/new_candidate3.apm'
V    = np.loadtxt(apm, usecols=(1,))
soun = np.loadtxt(apm, usecols=(0,), dtype=str)

##ICRF2 defining sources data file
icrf2f = '../data/icrf2-defining.dat'
defn = np.loadtxt(icrf2f, dtype=str, skiprows=20, usecols=(2,))

for i in range(len(soun)):
    if V[i] < 500:
        print(soun[i], file=fou)
    elif soun[i] in defn:
        print(soun[i], '  %8.2f'%V[i], file=fex)

fou.close()
fex.close()

print('Done!')