# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 15:20:29 2016

@author: neo
"""

#cat1 = 'new_candidate2.cat'
#cat1 = 'icrf1.cat'
cat1 = '394_sou.cat'
cat2 = 'icrf2.cat'
#cat3 = 'MFV247.cat'

f1 = open('../catalog/'+cat1,'r')
sou1 = f1.readlines()

f2 = open('../catalog/'+cat2,'r')
sou2 = f2.readlines()

#f3 = open('../catalog/'+cat3,'r')
#sou3 = f3.readlines()

f1.close()
f2.close()
#f3.close()

com = [ sou1[i] for i in range(len(sou1)) if sou1[i] in sou2 ]
#and sou1[i] in sou3]
print len(com)