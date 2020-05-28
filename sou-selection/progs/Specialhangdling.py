# -*- coding: utf-8 -*-
"""
Created on Tue May 24 14:20:14 2016

@author: Neo

Oct 25, 2016: updated by Niu.
"""

import numpy as np
import matplotlib.pyplot as plt

#ICRF2 defining sources data file
spef = '../list/39special.cat'
spen = np.loadtxt(spef, dtype=str)

V, V_err = np.loadtxt('../results/39special.apm', usecols=(1,4), unpack=True)

plt.plot(V_err, V, 'r+', markersize=15)
plt.xlabel('Linear drift $\mu (\mu as\,yr^{-1})$', fontsize = 15)
plt.ylabel('Uncertainty $\sigma _\mu(\mu as\,yr^{-1})$', fontsize = 15)
#plt.show()
plt.savefig('../plot/39Special.eps', dpi=100)
plt.close()