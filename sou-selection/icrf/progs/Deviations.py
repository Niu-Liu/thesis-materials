# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 20:43:50 2016

@author: neo

Calculate the standard variance, Allan Variance, goodness of fit.
"""

import fun
import numpy as np

def mean(x, sig):
    N = len(x)
    p = [sig[i]**-2 for i in range(N)]
    ave = sum([p[i]*x[i] for i in range(N)])/sum(p)
    sig_ave = np.sqrt(sum([p[i]**2*sig[i]**2 for i in range(N)])/\
    sum(p)**2)
#    x2 = [ x[i]**2 for i in range(N) ]
#    avex2 = sum([p[i]*x[i]**2 for i in range(N)])/sum(p)
#    dev = np.sqrt(avex2-ave**2)

#Another method to calculate Standard deviation
    dev = np.sqrt(sum([p[i]*(x[i]-ave)**2 for i in range(N)])/sum(p))
##
#These two methods equal in principle, but vary when calculating, may due to 
#the round off error of mechine.
    
    return [ave, sig_ave, dev]

##for test
#x = [1, 2]
#s = [0.4, 0.6]
#print mean(x,s)

#read the list of source name 
cat_dir = '../catalog/'
cat_fil = 'new_candidate1.cat'
fcat=open(cat_dir+cat_fil,'r')
catlg=fcat.readline()

dat_dir='../data/opa/'

out_fil = '../results/variance.dat'
fout = open(out_fil, 'w')

while len(catlg):
    soun=catlg.strip('\n')
    fdat=open(dat_dir+soun+'.dat','r')
    datcon=fdat.readlines()
    fdat.close()
    
    rac=[] #here it represent the curval segment
    dec=[]
    sigra=[]
    sigde=[]
    t=[]
    
    for j in range(3,len(datcon)-1):
        dat=datcon[j].split()
        alp = float(dat[1])
        det = float(dat[2])
        alpc = alp*np.cos(np.radians(det))
        
        t.append(fun.ADepo(float(dat[0])))     #yr
        rac.append(alpc*3.6e6)               #mas
        dec.append(det*3.6e6)               #mas
        sigra.append(float(dat[3])*np.cos(np.radians(det)))  #mas
        sigde.append(float(dat[4]))            #mas

# standard deviation        
    [ma, sigma, deva] = mean(rac, sigra) #mas
    [md, sigmd, devd] = mean(dec, sigde) #mas
    
#one-year-average
    ep1 = []
    ra1 = []
    sigra1 = []
    de1 = []
    sigde1 = []

    t1 = 1980
#    t1 = 1990
    t2 = 2016
    twd = t1
    sta = 0
    while twd < t2:
        ind = [ i for i in range(sta,len(t)) if twd<t[i]<twd+1 ]
        if len(ind):
            sta = ind[-1]+1
            [ma1, sigma1, deva1] = \
            mean(rac[ind[0]:sta], sigra[ind[0]:sta]) #mas
            [md1, sigmd1, devd1] = \
            mean(dec[ind[0]:sta], sigde[ind[0]:sta]) #mas
            
            ep1.append(twd)        #year
            ra1.append(ma1)        #mas
            sigra1.append(sigma1)  #mas
            de1.append(md1)        #mas
            sigde1.append(sigmd1)  #mas
            
        twd += 1
    
#Allan standard deviation
    N_ave = len(ep1)
    Ndat = 0
    AVar_a = 0.0
    AVar_d = 0.0
#WADEV, proposed by Malkin.    
    WAd_a = 0.0
    WAd_d = 0.0
    p_a = 0.0
    p_d = 0.0
    
    leap = []
    
    for i in range(len(ep1)-1):
        if ep1[i+1] - ep1[i] == 1:
            AVar_a += (ra1[i+1]-ra1[i])**2
            AVar_d += (de1[i+1]-de1[i])**2
            Ndat += 1
            
        pi_a = 1.0/(sigra1[i]**2 + sigra1[i+1]**2)
        pi_d = 1.0/(sigde1[i]**2 + sigde1[i+1]**2)
        WAd_a += pi_a*(ra1[i+1]-ra1[i])**2
        WAd_d += pi_d*(de1[i+1]-de1[i])**2
        
        p_a += pi_a
        p_d += pi_d
        
        leap.append(ep1[i+1] - ep1[i])
            
    AVar_a /= 2*Ndat
    AVar_d /= 2*Ndat
    
    ASd_a = np.sqrt(AVar_a)
    ASd_d = np.sqrt(AVar_d)
        
    WAd_a = np.sqrt(WAd_a/p_a/2.0)
    WAd_d = np.sqrt(WAd_d/p_d/2.0)
    
#store the eatimators into file.
#format:
#soure  Num_ave  Num_A           R.A.                 Decl.            Ob_gap
#                         mean  Std  AStd WAd | mean  Std  AStd WAd
#Num_ave: Number of one-year average points
#Num_A:   Number of data points to calculate Allan standard deviation 
#Ob_gap: Observation gap
#unit: mas for standard deviations
    print>>fout, soun + \
    " %2d %2d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %2d"\
    %(N_ave, Ndat, ma/3.6e6, deva, ASd_a, WAd_a, \
    md/3.6e6, devd, ASd_d, WAd_d, max(leap))
    catlg=fcat.readline()

fout.close()

print 'Done!'