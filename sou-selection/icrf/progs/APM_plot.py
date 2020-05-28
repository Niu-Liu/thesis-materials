# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 10:45:35 2016

@author: Neo
"""

import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Apm_Load 

dat_fil = '394_sou.apm'
#dat_fil = 'icrf2_sou.apm'
[Sou, V, V_err, V_RA, V_Dec, EV_RA, EV_Dec, RA, Dec] = Apm_Load(dat_fil)

#ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

Vdef = []
#Vwei = []
Vder = []
RAd = []

Vex = []
Ver = []
RAe = []

for i in range(len(Sou)):
    if Sou[i] in defn:
        Vdef.append(V[i])
#        Vwei.append(1.0/V_err[i]**2)
        Vder.append(V_err[i])
        RAd.append(RA[i])
        
    else:
        Vex.append(V[i])
        Ver.append(V_err[i])
        RAe.append(RA[i])

# Calculate the weighted average
#V = np.mat(V)
#Vwei = 1.0 / np.square(np.mat(V_err))
#ave = np.average(V, weights = Vwei)
#aved = np.average(Vdef, weights=Vwei)
print len(Vdef)

# Plot the histogram of the total apparent proper motion.   
numbins = [0, 20, 40, 60, 80, 100, 200, 1000, 5000]
count1, edge1= np.histogram(V, numbins)
count2, edge2= np.histogram(Vdef, numbins)
x1 = np.arange(len(count1))
x2 = x1+0.45
plt.bar(x1, count1, width=0.45, color='b', label='718 Candidates')
plt.bar(x2, count2, width=0.45, color='r', label='ICRF2 295 Defining')

for i in range(len(x1)):
    plt.text(x1[i]+0.15, count1[i], str(count1[i]), fontsize = 20)
    plt.text(x2[i]+0.15, count2[i], str(count2[i]), fontsize = 20)

yti = ['0', '20', '40', '60', '80', '100', '200', '1000']
plt.xticks(x1[:8], yti)
plt.xlabel('$\mu_{tol}(\mu as/yr$)', fontsize = 20)
plt.ylabel('Count', fontsize=15)
plt.legend(loc='upper right')
plt.show()
#plt.close()

#plot the total apparent proper motion VS RA
# This plot seems meaningless. -.-||
#plt.errorbar(RAd, Vdef, yerr = Vder, fmt='b.')
#plt.errorbar(RAe, Ver,  yerr = Ver,  fmt='r.')
#plt.xlim(0,360)
#plt.ylim(0,1000)
#plt.show()
#plt.close() 

# Plot Sources and its apparent proper motions.-
#flag1 = 0
#flag2 = 0
#for i in range(len(Sou)):
#    plt.arrow(RA[i], Dec[i], V_RA[i]/10, V_Dec[i]/10, width = 0.05, color = 'k')
#    if Sou[i] in defn:
#        plt.plot(RA[i], Dec[i], 'ro')
#    else:
#        plt.plot(RA[i], Dec[i], 'b.')
#
#plt.xlim(0, 360)
#plt.ylim(-90, 90)
#plt.xticks([0, 45, 90, 135, 180, 225, 270, 315, 360],\
#['0', '', '90', '', '180', '', '270', '', '360'])
#plt.yticks([-90, -60, -30, 0, 30, 60, 90])
#
#plt.xlabel('Right Ascension($\degree$)', fontsize = 20)
#plt.ylabel('  Declination  ($\degree$)', fontsize = 20)
#plt.grid(True)
#
### Proper motion scale
#plt.plot([10, 10], [-80.5, -79.5], 'k')
#plt.plot([20, 20], [-80.5, -79.5], 'k')
#plt.plot([10, 20], [-80  , -80  ], 'k')
#plt.text(8, -78, '$100\mu as/yr$', fontsize = 15)
#
#plt.legend()
#plt.show()
#plt.close()

# print the sources with APM >= 1000 microarcsecond/year.
#cat_dir = '../catalog/'
#cat_fil = '1mas.list'
#fout = file(cat_dir + cat_fil, 'w')
#for i in range(len(Sou)):
#    if V[i] >= 1000:
#        if Sou[i] in defn:
#            print >> fout, Sou[i]+'d'        
#        else:
#            print >> fout, Sou[i]
#fout.close()
