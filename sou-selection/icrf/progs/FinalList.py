# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 21:37:14 2016

@author: neo
"""

from Data_Load import Rank_Load 
[Sou, V, V_err, V_RA,V_Dec, Va_E, Vd_E, RA, Dec] = Rank_Load('GR.rank')
num = 394

fcat = open('../catalog/'+str(num)+'_sou.cat','w')
for i in range(num):
    print >>fcat, Sou[i]
fcat.close()