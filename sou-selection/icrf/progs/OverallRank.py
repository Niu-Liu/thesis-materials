# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:57:09 2016

To exclude the sources with apparent proper motion larger than 100 uas/yr
and do a simple rank

@author: Neo
"""

import numpy as np
from Data_Load import Apm_Load 

#dat_fil = '499_sou.apm'
dat_fil = '549_sou.apm'
[Sou, V, V_err, V_RA, V_Dec, Va_E, Vd_E, RA, Dec] = Apm_Load(dat_fil)

#cat_dir = '../catalog/'
#cat = 'NC3.cat'
#
#fcat = open(cat_dir + cat, 'w')

#dat = []
#for i in range(len(Sou)):
#    if V[i] < 1000:
#        dat.append((Sou[i], V[i], V_err[i],\
#        RA[i], V_RA[i], Va_E[i], Dec[i], V_Dec[i],Vd_E[i]))
#        print >> fcat, Sou[i]
#fcat.close()

dat = []
for i in range(len(Sou)):
    dat.append((Sou[i], V[i], V_err[i],\
        RA[i], V_RA[i], Va_E[i], Dec[i], V_Dec[i],Vd_E[i], V[i]/V_err[i]))
#fcat.close()

res_dir = '../results/'
res = 'OV.rank'
fres = open(res_dir + res,'w')

dtype = [('sou', 'S10'), ('u', float), ('u_e', float),\
        ('alp', float), ('ua', float), ('ua_e', float),\
        ('det', float), ('ud', float), ('ud_e', float), ('nu', float)]

a = np.array(dat, dtype=dtype)        
b = np.sort(a, order=['nu'])   

print len(b)    
       
for i in range(len(b)):
    print >> fres, b[i][0]+\
"   %10.4f    %10.4f  %10.4f    %10.4f  %10.4f  %10.4f  %10.4f  %10.4f"\
%(b[i][1],b[i][2],b[i][3],b[i][4],b[i][5],b[i][6],b[i][7],b[i][8])

fres.close()
print 'Done!'