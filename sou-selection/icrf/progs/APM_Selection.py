# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 21:58:57 2016

@author: neo
"""

from Data_Load import Apm_Load 
[Sou, V, V_err, V_RA, V_Dec, EV_RA, EV_Dec, RA, Dec] = Apm_Load('584_sou.apm')

##ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

m = 0
fcat = open('../catalog/new_candidate3.cat','w')
for i in range(len(Sou)):
    if V[i] < 500:
        print >>fcat, Sou[i]
        m += 1
    elif Sou[i] in defn:
        print Sou[i]
fcat.close()