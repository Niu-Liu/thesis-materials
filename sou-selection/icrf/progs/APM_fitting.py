# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 22:23:38 2016

Using the least square fitting to fit the apparent proper motion of Radio Sources

@author: Neo
"""

##Load the data
cat_dir = '../catalog/'
dat_dir = '../data/opa/'
out_dir = '../results/'

#sou_list = 'new_candidate1.cat' 
#out_fil = '613_sou.apm'

#sou_list = 'new_candidate2.cat' 
#out_fil = '584_sou.apm'

#sou_list = 'new_candidate3.cat' 
#out_fil = '549_sou.apm'

#sou_list = '322_sou.cat' 
#out_fil = '322_sou.apm'

#sou_list = '351_sou.cat' 
#out_fil = '351_sou.apm'

sou_list = '394_sou.cat' 
out_fil = '394_sou.apm'

##generate the apparent proper motion data of some specific catalogues.
##for ICRF1
#out_fil = 'icrf1_sou.apm'
#sou_list = 'icrf1.cat'
##for ICRF2
#out_fil = 'icrf2_sou.apm'
#sou_list = 'icrf2.cat'
#for MFV247
#out_fil = 'MFV247_sou.apm'
#sou_list = 'MFV247.cat'

#out_fil = 'AMS260_sou.apm'
#sou_list = 'AMS260.cat'

import fun
import matrix_weighed
from numpy import *
from math import *

# read catalog file to get name of sources.
fcat = open(cat_dir + sou_list,'r')
catcon = []
st = fcat.readline()
while len(st):
    catcon.append(st)
    st = fcat.readline()
fcat.close()

fout = open(out_dir + out_fil, 'w')

for i in range(len(catcon)):
    sou_name = str(catcon[i]).strip('\n')
    fdat = open(dat_dir + sou_name + '.dat','r')
    filcon = fdat.readlines()
    datcon = filcon[3:]
    N = len(datcon)
    
    t=[]
    RA=[]
    DEC=[]
    ERA=[]
    EDEC=[]

#data format
#Obs_moment   R.A.(o)   Dec.(o)  Err. (mas)  Corr. Del.  IERS   IVS   
#Unit: year for time and microarcsecond for angles 
    for j in range(N):
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
    
#    r1 = RA  - t*mu_ra - ra0
#    r2 = DEC - t*mu_ra - de0
#    e1 = [sqrt((emu_ra*t[i])**2 + era0**2) for i in range(N)]
#    e2 = [sqrt((emu_ra*t[i])**2 + era0**2) for i in range(N)]
    
    mu_ra = mu_ra*cos(radians(de0))
    emu_ra = emu_ra*cos(radians(de0))
    
    mu_tol = sqrt(mu_ra**2 + mu_de**2)
    emu_tol = sqrt(emu_ra**2 + emu_de**2)
# Write the results into outpu file
##Unit: microarcsecond/year for proper motion, and degree for angles 
#data format
#source   P.M.(uas/yr)  Err. (uas/yr)   R.A.(o)   Dec.(o)  Err. (mas)
#            (Total R.A. Dec.)
#R.A. and Dec. corrspond to JD2000.0
    print >> fout, sou_name + "  %10.4f    %10.4f  %10.4f    %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f" \
    %(mu_tol,mu_ra, mu_de, emu_tol, emu_ra, emu_de,\
    ra0, de0, era0, ede0)
    
fout.close()
print 'Done!'                                           
