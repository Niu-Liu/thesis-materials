# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 16:19:11 2016

@author: Neo
"""

#import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Apm_Load 

##ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

dat_fil = '499_sou.apm'
[Sou, V, V_err, V_RA, V_Dec, RA, Dec] = Apm_Load(dat_fil)

Vdef = []
#Vwei = []
Vder = []
RAd = []
Decd = []

Vex = []
Ver = []
RAe = []
Dece = []

for i in range(len(Sou)):
    if Sou[i] in defn:
        Vdef.append(V[i])
#        Vwei.append(1.0/V_err[i]**2)
        Vder.append(V_err[i])
        RAd.append(RA[i])
        Decd.append(Dec[i])
    else:
        Vex.append(V[i])
        Ver.append(V_err[i])
        RAe.append(RA[i])
        Dece.append(Dec[i])
        
#plt.plot(Vder, Vdef, 'r.')
#plt.plot(Ver, Vex, 'b.')
##plt.ylim([0, 1000])
##plt.xlim([0,  200])
#plt.show()

#Normalized Linear drift
Nld = [ V[i]/V_err[i] for i in range(len(Sou)) ]
#plt.plot(V, Nld, '.')
plt.plot(Dec, Nld, '.')
#plt.xlim([0, 100])
plt.show()