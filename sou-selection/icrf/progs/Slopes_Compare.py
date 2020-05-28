#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 19:03:15 2015

@author: Neo

Compare the slopes of coordinates datas smoothed using different ways with each
 others.

"""

import math

## Fetch apparent proper motion data.
list_dir='../results/'
listf='sou.apm.Nov04'
fdat=open(list_dir+listf,'r')
dat=fdat.readlines()
dat=dat[2:len(dat)]
fdat.close()

#[sou,u,ua,ud,sigua,sigud,alp0,det0,siga0,sigd0,chia,chid]
sou=[]
ua=[]
ud=[]
alp0=[]
det0=[]

for i in range(len(dat)):
    datlin=dat[i].split()
    sou.append(datlin[0])
    ua.append(float(datlin[2]))
    ud.append(float(datlin[3]))
    alp0.append(float(datlin[7]))
    det0.append(float(datlin[8]))
    
list1='sou.apm1990'
fdat=open(list_dir+list1,'r')
dat=fdat.readlines()
dat=dat[2:len(dat)]
fdat.close()

out_dir='../results/'
outf='apm.dif'+list1[7:len(list1)]
fp=open(out_dir+outf,'w')
print >>fp, "#sources    u_alp*cos(det)    u_det   alp "

for i in range(len(dat)):
    datlin=dat[i].split()
    p=sou.index(datlin[0])
    det1=-ua[p]+float(datlin[2])
    det2=-ud[p]+float(datlin[3])
    print >>fp,  sou[p]+"  %10.4f    %10.4f  %10.4f  %10.4f  %6.2f" \
    %(det1*math.cos(det0[p]/180.0*math.pi), det1/math.fabs(ua[p])*100,\
    det2, det2/math.fabs(ud[p])*100, alp0[p])

fp.close()

print 'Done!'
