# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 09:16:33 2016

Calculation the Allan variance
!!!!!
This code seems to have some issues.
@author: Neo
"""

dat_dir = '../data/Annual_Smooth/'
cat_dir = '../catalog/'
res_dir = '../results/'

cat = 'new_candidate2.cat'
res = '512Allan.var'

fcat = open(cat_dir + cat,'r')
fres = open(res_dir + res,'w')
soun = fcat.readline().strip('\n')

N = 36

while len(soun):
    fdat = open(dat_dir + soun + '.dat','r')
    dat  = fdat.readline().split()
    ra  = []
    dec = []
    
    while len(dat):
        ra.append(float(dat[0]))
        dec.append(float(dat[2]))
        dat = fdat.readline().split()
        
    if len(ra) == N:
        x1 = 0
        x2 = 0
        for i in range(N-1):
            x1 += (ra[i+1] - ra[i])**2
            x2 += (dec[i+1] - dec[i])**2
            
        x1 = (x1/(2.0*N))**0.5
#        /3.6e6
        x2 = (x2/(2.0*N))**0.5
#        /3.6e6
        
        print>>fres, soun + "    %6.4f    %6.4f"%(x1,x2)        

    else:
        print soun + "ERROR!"
    
        
    soun = fcat.readline().strip('\n')
#    soun =[]
    fdat.close()    
    
fres.close()
fcat.close()

print 'Done!'
