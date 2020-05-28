# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 09:35:15 2016

@author: Neo
"""

cat_dir = '../catalog/'
#for icrf1
inp = 'icrf1.def'
oup = 'icrf1.cat'
#for icrf2
inp = 'icrf2-defining.dat'
oup = 'icrf2.cat'

fin = open(cat_dir + inp, 'r')
fou = open(cat_dir + oup, 'w')

lin =fin.readline()
while len(lin):
    if lin[0:4] == 'ICRF':
#        for icrf1
#        print>>fou, lin[24:32]
#        for icrf2
        print>>fou, lin[23:31]
    lin =fin.readline()
    
fin.close()
fou.close()
print 'Done!'
