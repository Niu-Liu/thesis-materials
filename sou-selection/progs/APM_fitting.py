# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 22:23:38 2016

Using the least square fitting to fit the apparent proper motion of Radio Sources

@author: Neo

Oct 17, 2016: updated by Niu.
"""

import numpy as np
from LinearFitting import APM
import time
from fun import ADepoA
cos = np.cos

##Load the data
cat_dir = '../list/'
dat_dir = '../data/opa/'
out_dir = '../results/'

#sou_list = 'opa.list' 
#out_fil = 'opa_all.apm'
#
#sou_list = 'new_candidate4.cat' 
#out_fil  = 'new_candidate4.apm'

sou_list = '39special.cat' 
out_fil  = '39special.apm'

## read catalog file to get name of sources.
soun = np.loadtxt(cat_dir + sou_list, dtype=str)
## output file. 
fout = open(out_dir + out_fil, 'w')
print(time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\
    time.localtime(time.time())), file=fout)
print('## source list: %s'%sou_list, file=fout)
print('''##data format
##source   P.M.(uas/yr)  Err. (uas/yr)   R.A.(o)   Dec.(o)  Err. (mas)
##            (Total R.A.* Dec.)
##R.A. and Dec. corrsponding to JD2000.0
############################################################################''', file=fout)

for i in range(len(soun)):
    sou_name = soun[i]
#    print sou_name 
##data format
##epo(Julian day)   R.A.(o)   Dec.(o)  Err. (mas)  Corr. Del.  IERS   IVS   
    epo, ra, dec, era, edec = \
    np.loadtxt(dat_dir+sou_name+'.dat', usecols=list(range(5)), unpack=True)
    ## Unit: Julian day, deg, deg, mas, mas
    
    if epo.size>1:
        epo = ADepoA(epo)                       #yr
        ra  = ra*3.6e6                          #mas
        dec = dec*3.6e6                         #mas 

        [pmRA, pmDE, RA0, DE0, epmRA, epmDE, eRA0, eDE0] = \
        APM(epo, ra, era, dec, edec)
        
        era0, ede0 = eRA0, eDE0
##  mas -> deg    
        ra0, de0 = RA0/3.6e6, DE0/3.6e6
##  mas/yr -> uas/yr
        mu_ra, mu_de = pmRA*1.0e3*cos(np.deg2rad(de0)), pmDE*1.0e3
        emu_ra, emu_de = epmRA*1.0e3*cos(np.deg2rad(de0)), epmDE*1.0e3
    
    else:
        mu_ra, mu_de = 0.0, 0.0
        emu_ra, emu_de = 1.0e6, 1.0e6
        ra0, era0 = ra,  era 
        de0, ede0 = dec, edec
        
    mu_tol  = np.sqrt( mu_ra**2 +  mu_de**2)
    emu_tol = np.sqrt(emu_ra**2 + emu_de**2)
    
    earray = np.array([emu_ra, emu_de, emu_tol])
    emu_ra, emu_de, emu_tol = np.where(earray>1.0e6, 1.0e6, earray)
# Write the results into outpu file
##Unit: microarcsecond/year for proper motion, and degree for angles 
#data format
#source   P.M.(uas/yr)  Err. (uas/yr)   R.A.(o)   Dec.(o)  Err. (mas)
#            (Total R.A. Dec.)
#R.A. and Dec. corrspond to JD2000.0
    print(sou_name + ("  %12.4f"*6+"  %14.10f"*2+"  %12.4f"*2) \
    %(mu_tol,mu_ra, mu_de, emu_tol, emu_ra, emu_de,\
    ra0, de0, era0, ede0), file=fout)
    
fout.close()
print('Done!')                                           