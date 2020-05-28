# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:52:15 2015

@author: Neo
"""

#select_on_pos.py

import fun
#import numpy as np
datdir='../results/'
#datf='sou.apm.Nov04SAL'
#datf='sou.apm.Nov04SA'
datf='sou.apm1990'

#outdir='../catalog/'
#outf='SAL.list'

## Fetch apparent proper motion data.
fdat=open(datdir+datf,'r')
dat=fdat.readlines()
dat=dat[2:len(dat)]
fdat.close()

## ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2,'r')
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

#[sou,u,ua,ud,sigua,sigud,alp0,det0,siga0,sigd0,chia,chid]
L=len(dat)
sou=[]
u=[]
sigu=[]
ua=[]
ud=[]
sigua=[]
sigud=[]
alp0=[]
det0=[]

for i in range(len(dat)):
    datlin=dat[i].split()
    sou.append(datlin[0])
    u.append(float(datlin[1]))
    ua.append(float(datlin[2]))
    ud.append(float(datlin[3]))
    sigu.append(float(datlin[4]))
    sigua.append(float(datlin[5]))
    sigud.append(float(datlin[6]))
    alp0.append(float(datlin[7]))
    det0.append(float(datlin[8]))
    
# nodes which dic=vide radio sources into 6 groups.    
#node=fun.sou_division(det0,6)
#
#print ' the division nodes of sources are :', node

# get the mean apparent proper motion volicity of ICRF2.
#meanu=fun.mean_apm(sou,defn,u,sigu)
#print meanu

# Calculate axis rotation rate of the ICRF2 295 radio sources.
i=0
while i< len(sou):
    if sou[i] in defn:
        i+=1
    else:
        del sou[i], ua[i], ud[i], sigua[i], sigud[i], alp0[i], det0[i]
x=fun.axis_rot(ua,ud,sigua,sigud,alp0,det0)
print x

## Calculate axis rotation rate of the whole 691 radio sources.
#x=fun.axis_rot(ua,ud,sigua,sigud,alp0,det0)
#print x

# get the mean apparent proper motion volicity of whole 695 radio sources.
#meanu=fun.mean_apm(sou,sou,u,sigu)
#print meanu

#j=0
#for i in range(len(sou)):
#    if u[i]<=21.97:
#        j+=1
#
#print j
