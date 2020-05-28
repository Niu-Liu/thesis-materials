# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:27:29 2016

@author: Neo
"""

import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Apm_Load 

dat_fil = '499_sou.apm'
[Sou, V, V_err, V_RA, V_Dec, RA, Dec] = Apm_Load(dat_fil)

## I divide the sources into 4 groups according to its declination.
## Declination: -90 ~ -30; -30 ~ 0; 0 ~ 30; 30 ~ 90
V1 = []
V2 = []
V3 = []
V4 = []

dat1 = []
dat2 = []
dat3 = []
dat4 = []


for i in range(len(Sou)):
    if Dec[i] < -30:
        V1.append(V[i])
        dat1.append([Sou[i], V[i]])
    elif Dec[i] < 0:
        V2.append(V[i])
        dat2.append([Sou[i], V[i]])
    elif Dec[i] < 30:
        V3.append(V[i])
        dat3.append([Sou[i], V[i]])
    else:
        V4.append(V[i])
        dat4.append([Sou[i], V[i]])
        
#Plot the counts of total APM.
numbins = [0, 20, 40, 60, 80, 100, 200, 1000, 5000]
count1, edge1= np.histogram(V1, numbins)
count2, edge2= np.histogram(V2, numbins)
count3, edge3= np.histogram(V3, numbins)
count4, edge4= np.histogram(V4, numbins)

x1 = np.arange(len(count1))
x2 = x1+0.24
x3 = x2+0.24
x4 = x3+0.24
plt.bar(x1, count1, width=0.24, color='b', label='$-90\degree \sim -30\degree$')
plt.bar(x2, count2, width=0.24, color='r', label='$-30\degree \sim   0\degree$')
plt.bar(x3, count3, width=0.24, color='k', label='$  0\degree \sim  30\degree$')
plt.bar(x4, count4, width=0.24, color='y', label='$ 30\degree \sim  90\degree$')

for i in range(len(x1)):
    plt.text(x1[i]+0.05, count1[i], str(count1[i]), fontsize = 20)
    plt.text(x2[i]+0.05, count2[i], str(count2[i]), fontsize = 20)
    plt.text(x3[i]+0.05, count3[i], str(count3[i]), fontsize = 20)
    plt.text(x4[i]+0.05, count4[i], str(count4[i]), fontsize = 20)

yti = ['0', '20', '40', '60', '80', '100', '200', '1000']
plt.xticks(x1[:8], yti)
plt.xlabel('$\mu_{tol}(\mu as/yr$)', fontsize = 20)
plt.ylabel('Count', fontsize=15)
plt.legend(loc='upper right')
plt.show()
#plt.close()
print len(V1), len(V2), len(V3), len(V4)

#Sort dat, and write a rank list 
dat1.sort(key=lambda x:x[1])
dat2.sort(key=lambda x:x[1])
dat3.sort(key=lambda x:x[1])
dat4.sort(key=lambda x:x[1])

#if we need 300 candidates, then 75 for each group.
dat1 = dat1[:75]
dat2 = dat2[:75]
dat3 = dat3[:75]
dat4 = dat4[:75]

out_dir = '../catalog/'
out_fil = '300can.list'
fout = file(out_dir + out_fil, 'w')
for i in range(len(dat1)):
    print >> fout, dat1[i][0]
for i in range(len(dat2)):
    print >> fout, dat2[i][0]
for i in range(len(dat3)):
    print >> fout, dat3[i][0]
for i in range(len(dat4)):
    print >> fout, dat4[i][0]
fout.close()
