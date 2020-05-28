# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 16:37:21 2016

@author: Neo

regeneration of candidates
"""

from Data_Load import Apm_Load 

dat_fil = '635_sou.apm'

[Sou, V, V_err, V_RA, V_Dec, RA, Dec] = Apm_Load(dat_fil)

cat_dir = '../catalog/'
can = 'cadidates2.cat'

sou = []
for i in range(len(Sou)):
    if V[i] < 100:
        sou.append(Sou[i])
        
print len(sou)
