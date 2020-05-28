# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 16:35:47 2016
7 defining sources is ejected by selection algorithm.
Here I just want to check their time series information.
@author: Neo
"""

ex = ['0458+138', '1448-648', '1554-643', '1611-710', \
      '1633-810', '1824-582', '2344-514']
      
dat_dir = '../data/opa/'

import fun
import matrix_weighed
from numpy import *
from math import *

for i in range(len(ex)):
    sou = ex[i]
    fdat = open(dat_dir + sou + '.dat','r')
    filcon = fdat.readlines()
    datcon = filcon[3:]
    
    t=[]
    RA=[]
    DEC=[]
    ERA=[]
    EDEC=[]

#data format
#Obs_moment   R.A.(o)   Dec.(o)  Err. (mas)  Corr. Del.  IERS   IVS   
#Unit: year for time and microarcsecond for angles 
    for j in range(len(datcon)):
        data = str(datcon[j]).split()[:5]
        t.append( fun.ADepo(float(data[0])) )
        RA.append( float(data[1]) * 3.6e9 )
        DEC.append( float(data[2]) * 3.6e9 )
        ERA.append( float(data[3]) * 1.0e3 )
        EDEC.append( float(data[4]) * 1.0e3 )
        
    t = array(t)
    RA = array(RA)
    DEC = array(DEC)
    ERA = array(ERA)
    EDEC = array(EDEC)
    
    [ra,sigra,de,sigde] = \
    matrix_weighed.matrix_weighed_fitting(RA,ERA,DEC,EDEC,t)
    
    [mu_ra,ra0] = [ra[0], ra[1]/3.6e9]
    [emu_ra,era0] = [sigra[0], sigra[1]/1.0e3]
    [mu_de,de0] = [de[0], de[1]/3.6e9]
    [emu_de,ede0] = [sigde[0], sigde[1]/1.0e3]
    
    mu_ra = mu_ra*cos(radians(de0))
    emu_ra = emu_ra*cos(radians(de0))
    
    mu_tol = sqrt(mu_ra**2 + mu_de**2)
    emu_tol = sqrt(emu_ra**2 + emu_de**2)

    print sou + "  %10.4f    %10.4f  %10.4f    %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f" \
    %(mu_tol,mu_ra, mu_de, emu_tol, emu_ra, emu_de,\
    ra0, de0, era0, ede0)
    
fout.close()
print 'Done!'          
