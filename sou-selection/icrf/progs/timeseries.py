#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 20:08:53 2015

@author: Neo
"""

#!/bin/python
#

# data and catalog file and path
#dat_dir='../data/opa/Raw/'
dat_dir='../data/opa/Post1990.0/'
cat_dir='../catalog/'
cat_fil='new_candidate1990.list'

#dat_dir='../data/Analysis/Annual_Smooth/'
#cat_dir='../catalog/'
#cat_fil='A_smth.list'
#
#dat_dir='../data/Analysis/Semi-Annual_Smooth/'
#cat_dir='../catalog/'
#cat_fil='SA_smth.list'

#dat_dir='../data/Analysis/Semi-Annual_Smooth_L/'
#cat_dir='../catalog/'
#cat_fil='SAL_smth.list'

#string array of candidate radio sources names
#fcat=open(cat_dir+cat_fil,'r')
#catlg=fcat.readline()
##catlg=fcat.readlines()
#fcat.close()
#catlg=catlg.split('\t')
#del catlg[len(catlg)-1]

fcat=open(cat_dir+cat_fil,'r')
catlg=fcat.readlines()
fcat.close()

#result file, apparent proper motion of radio sources
res_fil='sou.apm1990'
fres=open(res_fil,'w')
fres.write('##Sources \t mu_tol \t mu_RA \t mu_DE \t \tsigma \t RA \t DE \t '+\
'\tsigma \t chi-square\n')
fres.write('## \t\t uas/yr \t uas/yr \t \tuas/yr \t\t deg \t deg \t'+\
'\tdeg \n')

import fun
t0=fun.ADepo(51544.5)
#import numpy as np

for i in range(len(catlg)):
    sou=catlg[i].strip('\n')
    
#    fdat=open(dat_dir+sou+'.txt','r')       #1
#    datcon=fdat.readlines()
#    leng=len(datcon)-3                      #1
    
    fdat=open(dat_dir+sou+'.dat','r')      #2
    datcon=fdat.readlines()

    
    ra=[]
    de=[]
    sigra=[]
    sigde=[]
    t=[]
#    
#    for j in range(3,len(datcon)-1):
##    for j in range(len(datcon)):
#        dat=datcon[j].split()
#        t.append(fun.ADepo(float(dat[0]))-t0)  #yr
##        t.append(float(dat[0])-t0)  #yr
#        ra.append(float(dat[1])*3.6e9)         #uas
#        de.append(float(dat[2])*3.6e9)         #uas
#        sigra.append(float(dat[3])*1e3)        #uas
#        sigde.append(float(dat[4])*1e3)        #uas
#    
    for j in range(len(datcon)):
        dat=datcon[j].split()
        t.append(float(dat[0])-t0)             #yr
        ra.append(float(dat[1]))         #uas
        de.append(float(dat[3]))         #uas
        sigra.append(float(dat[2]))        #uas
        sigde.append(float(dat[4]))        #uas
##        
        
    [u,sigu,ua,ud,sigua,sigud,alp0,det0,siga0,sigd0,chia,chid]=\
    fun.pm(ra,sigra,de,sigde,t,1979.0-t0)
    
#    fres.write(sou+'\t'+str(u)+'\t'+str(ua)+'\t'+str(ud)+'\t'+str(sigua)+'\t'+str(sigud)+\
#    '\t'+str(alp0)+'\t'+str(det0)+'\t'+str(siga0)+'\t'+str(sigd0)+'\t'+\
#    str(chia)+'\t'+str(chid)+'\n')
    print >> fres, sou+"  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %6.2f  %6.2f  %6.2f  %6.2f  %10.4f  %10.4f" \
    %(u,ua,ud,sigu,sigua,sigud,alp0,det0,siga0,sigd0,chia,chid)
    
fres.close() 

print 'Done!'
