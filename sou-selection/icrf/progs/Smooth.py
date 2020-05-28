#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
smooth.py

Created on Sun Nov  1 22:41:31 2015

@author: Neo
"""

import math
import numpy as np
import fun

def Duplicate_removal(ids):
    func = lambda x,y:x if y in x else x + [y]
    return reduce(func, [[], ] + ids)

def mavr(x,sig):
#   calculate the weighed mean and its uncertainty.
    p=[ sig[i]**-2 for i in range(len(sig)) ]
    px=[ p[i]*x[i] for i in range(len(p)) ]
    wrm= 1.0/sum(p)
    mean= sum(px)*wrm
    return [mean,math.sqrt(wrm)]
    
#def annual_smooth(alp,siga,det,sigd,epo,soun):
##way 1:   Annual smmoth of data points.
#    dest_dir='../data/Annual_Smooth/'
##    dest_dir='../data/Analysis/Semi-Annual_Smooth/'
##    dest_dir='../data/Analysis/Semi-Annual_Smooth_L/'
#    smth_fil=soun+'.dat'
#    fd=open(dest_dir+smth_fil,'w')
#    
##   if the number of sessions within a windows is less than 3,    
##   set data point to data0 with high uncertainty sigma0
##   data0=0.0
##   sigma0=1e10    
##   But i'd like to set it as empty point.    
#    
##   epoch span, length of window and step.
#    epo_b = 1979.0
#    epo_e = 2016.0
##   Annual Smooth
#    widl=1.0
#    step=1.0
#    
##Way 2: Semi_annual Smooth
##    widl=0.5
##    step=0.5  
#    
##Way 3: Semi_annual Smooth, Lambert & Gontier, A&A, 2009
##    widl = 2
##    step = 0.5  
#    
#    wid = epo_b
##    while wid < epo_e:
##        dat=[ [alp[i],siga[i],det[i],sigd[i],epo[i]] for i in range(len(epo)) \
##        if wid<=epo[i]<=wid+widl ]
##        if len(dat)>=3:
##            dat=np.matrix(dat)
##            [alp1,siga1]=mavr(dat[:,0],dat[:,1])
##            [det1,sigd1]=mavr(dat[:,2],dat[:,3])
##            print >> fd, "%f %.3f %.1f %.3f %.1f" % (sum(dat[:,4])/len(dat),\
##            alp1,siga1,det1,sigd1)
##            
##        wid+=step
#        
##just for way3
#    alp0=[]
#    siga0=[]
#    det0=[]
#    sigd0=[]
#    t0=[]
#    
#    while wid < epo_e:
#        dat=[ [alp[i],siga[i],det[i],sigd[i],epo[i]] for i in range(len(epo)) \
#        if wid <= epo[i] <= wid+widl ]
#        if len(dat) >= 3:
#            dat=np.matrix(dat)
#            [alp1,siga1]=mavr(dat[:,0],dat[:,1])
#            [det1,sigd1]=mavr(dat[:,2],dat[:,3])
#            t0.append(sum(dat[:,4])/len(dat))
#            alp0.append(alp1)
#            siga0.append(siga1)
#            det0.append(det1)
#            sigd0.append(sigd1)
#            
#        else:       
#            t0.append(wid)
#            alp0.append(0.0)
#            siga0.append(1.0e10)
#            det0.append(0.0)
#            sigd0.append(1.0e10)
#            
#            
#        wid+=step
#        
##    t0=Duplicate_removal(t0)
##    alp0=Duplicate_removal(alp0)
##    siga0=Duplicate_removal(siga0)
##    det0=Duplicate_removal(det0)
##    sigd0=Duplicate_removal(sigd0)
#        
#    i=0
#    while i<len(t0)-1:
#        if t0[i] == t0[i+1]:
#            del t0[i+1],alp0[i+1],siga0[i+1],det0[i+1],sigd0[i+1]
#        else:
#            print >> fd, "%f %.3f %.1f %.3f %.1f" % (t0[i],\
#    alp0[i],siga0[i],det0[i],sigd0[i])
#            i+=1
#        
#    fd.close()
    
    
#rewritten by N. Liu
#Thu Mar 10 20:58:38 CST 2016
def annual_smooth(alp,siga,det,sigd,epo,soun,fin,fex):    
#   if the number of sessions within a windows is less than 3,    
#   set data point to data0 with high uncertainty sigma0
#   data0=0.0
#   sigma0=1e10  
#   Here we set a flag "Flag" to record the length of successive years 
#   within  number of observation sessions is less 3
#   if Flag >= 3
#   then the sources is not suitable because of bad continuity

    Flag  = 0
    count = 0
#   epoch span, length of window and step.
    epo_b = 1980.0
    epo_e = 2016.0

    widl=1.0
    step=1.0
    
    wid = epo_b
    alp0=[]
    siga0=[]
    det0=[]
    sigd0=[]
    t0=[]
    
    tag = 1    
    
    while wid < epo_e:
        dat=[ [alp[i],siga[i],det[i],sigd[i],epo[i]] for i in range(len(epo)) \
        if wid <= epo[i] <= wid+widl ]
        if len(dat) >= 3:
            dat=np.matrix(dat)
            [alp1,siga1]=mavr(dat[:,0],dat[:,1])
            [det1,sigd1]=mavr(dat[:,2],dat[:,3])
            t0.append(wid)
            alp0.append(alp1)
            siga0.append(siga1)
            det0.append(det1)
            sigd0.append(sigd1)
            
            Flag  = 1
            count = 0
            
        else:       
            t0.append(wid)
            alp0.append(0.0)
            siga0.append(1.0e10)
            det0.append(0.0)
            sigd0.append(1.0e10)
            
            if Flag:
                count += 1
#            count += 1
        
        if count > 10:
            tag = 0
               
        wid += step
        
    
    dest_dir='../data/Annual_Smooth/'
    smth_fil=soun+'.dat'
    fd=open(dest_dir+smth_fil,'w')
        
    for i in range(len(t0)):
#            print >> fd, "%f %.3f %.1f %.3f %.1f" % (t0[i],\
#                    alp0[i],siga0[i],det0[i],sigd0[i])
#There is no need to write the time information into the datafile
        print >> fd, "%.6f %.4f %.6f %.4f" % (\
    alp0[i],siga0[i],det0[i],sigd0[i])
                    
    fd.close()
        
    if tag:
        print >> fin, soun
        
    else:
        print >> fex, soun
        
        if soun in defn:
            print soun
        
##ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()    
        
# path of folders and files
dat_dir='../data/opa/'
cat_dir='../catalog/'
cat_fil='new_candidate1.cat'

#string array of candidate radio sources names
fcat=open(cat_dir+cat_fil,'r')
catlg=fcat.readline()

new_cat = 'new_candidate2.cat'
ex_cat  = 'contunuity_ex.cat'

fnew = open(cat_dir + new_cat, 'w')
fex  = open(cat_dir + ex_cat,  'w')

#for i in range(len(catlg)):
while len(catlg):
    soun=catlg.strip('\n')
    fdat=open(dat_dir+soun+'.dat','r')
    datcon=fdat.readlines()
    
    rac=[] #here it represent the curval segment
    de=[]
    sigra=[]
    sigde=[]
    t=[]
    
    for j in range(3,len(datcon)-1):
        dat=datcon[j].split()
        alp = float(dat[1])
        det = float(dat[2])
        alpc = alp*math.cos(math.radians(det))
        
        t.append(fun.ADepo(float(dat[0])))     #yr
        rac.append(alpc*3.6e6)               #deg
        de.append(det*3.6e6)               #deg
        sigra.append(float(dat[3]))            #uas
        sigde.append(float(dat[4]))            #uas
        
    annual_smooth(rac,sigra,de,sigde,t,soun,fnew,fex)
    
    catlg=fcat.readline()
 
fcat.close()
fnew.close()
fex.close()
   
print 'Done!'
