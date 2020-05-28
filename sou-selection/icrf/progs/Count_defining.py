# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:32:24 2016

@author: neo

count the number of defining sources in a given list
"""

##ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

#dat = '../catalog/new_candidate2.cat'
dat = '../catalog/contunuity_ex.cat'

fdat = open(dat, 'r')
st = fdat.readline().strip('\n')
cou = 0
N = len(st)
while len(st):
    if st in defn:
        cou += 1
        
    st = fdat.readline().strip('\n') 
    
print 'defining sources    '
print cou
