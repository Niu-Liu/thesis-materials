# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 10:28:00 2016

Get the parallax information from TGAS data for O-B5 stars in Tycho-2

@author: Neo
"""

import numpy as np
import time

dat_dir = '/Users/Neo/Astron/Work/201609_HCRF/data/'

fty2 = dat_dir + 'tycho2_III231.mk'
ftgast= dat_dir + 'TagsSource.tyc'
ftgash= dat_dir + 'TagsSource.hip'
tyc2t = dat_dir + 'tyc2.tyc'
tyc2h = dat_dir + 'tyc2.hip'
outh  = dat_dir + 'Tags_hip_III231_mk.dat'
outt  = dat_dir + 'Tags_tyc_III231_mk.dat'

## Mapping between of Hipparcos and Tycho-2 identifiers. 
Hip, TYC = np.loadtxt(dat_dir+'Hip_Tyc.map', dtype=str, unpack=True)
    
flog = open('../logs/cross_ident.log', 'a')
print>>flog,'#################################################################'
print>>flog, ''' Cross identification between TGAS and III231.ob begins.'''
print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    
taggt = np.loadtxt(ftgast, dtype=str, usecols=(0,))
taggh = np.loadtxt(ftgash, dtype=str, usecols=(0,))
   
tagtt = np.loadtxt(tyc2t, dtype=str, usecols=(0,))
tagth = np.loadtxt(tyc2h, dtype=str, usecols=(0,))
    
## data from Tycho-2
rath, decth, pmrath, pmdecth, ra_errth, dec_errth, pmra_errth, pmdec_errth = \
    np.loadtxt(tyc2h, usecols=(range(1,9)), unpack=True)
ratt, dectt, pmratt, pmdectt, ra_errtt, dec_errtt, pmra_errtt, pmdec_errtt = \
    np.loadtxt(tyc2t, usecols=(range(1,9)), unpack=True)
    
## data from TGAS 
ragh, ra_errgh, decgh, dec_errgh, plxh, plx_errh, \
    pmragh, pmra_errgh, pmdecgh, pmdec_errgh, lh, bh \
    = np.loadtxt(ftgast, usecols=(range(1,11)+[21, 22]), unpack=True)    
ragt, ra_errgt, decgt, dec_errgt, plxt, plx_errt, \
    pmragt, pmra_errgt, pmdecgt, pmdec_errgt, lt, bt  \
    = np.loadtxt(ftgast, usecols=(range(1,11)+[21, 22]), unpack=True)


num = 0     
N = 0
print 'cross idenfication begins! '
    
fouh = open(outh, 'w')
print>>fouh, '''## Information for O-B5 stars of Hipparcos 
##Data souces: position and pm from tycho-2; plx data from TGAS
##Hip_identifier  ra(deg)   e_ra(mas)   dec(deg)   e_dec(mas)   pmra(mas/yr)   pmdec(mas/yr)   e_pmra(mas/yr)   e_pmdec(mas/yr)  Parallax(mas)   e_plx(mas)   VTmag  BTmag  Teff(K)  SpClass  LClass'''

'''## Information for O-B5 stars of Hipparcos 
##Data souces: T for Tycho-2; G for TGAS; plx data from TGAS
##Hip_identifier  raT(deg)   e_raT(mas)   decT(deg)   e_decT(mas)   pmraT(mas/yr)   pmdecT(mas/yr)   e_pmraT(mas/yr)   e_pmdecT(mas/yr)  raG(deg)   e_raG(mas)   decG(deg)   e_decG(mas)   pmraG(mas/yr)   pmdecG(mas/yr)   e_pmraG(mas/yr)   e_pmdecG(mas/yr)   Parallax(mas)   e_plx(mas)   VTmag  BTmag  Teff(K)  SpClass  LClass'''

fout = open(outt, 'w')
print>>fout, '''## Information for O-B5 stars of Hipparcos 
##Data souces: position and pm from tycho-2; plx data from TGAS
##Tycho-2_identifier  ra(deg)   e_ra(mas)   dec(deg)   e_dec(mas)   pmra(mas/yr)   pmdec(mas/yr)   e_pmra(mas/yr)   e_pmdec(mas/yr)  Parallax(mas)   e_plx(mas)   VTmag  BTmag  Teff(K)  SpClass  LClass'''


'''## Information for O-B5 stars of Tycho-2 
##Data souces: T for Tycho-2; G for TGAS; plx data from TGAS
## Tycho-2_identifier  raT(deg)   e_raT(mas)   decT(deg)   e_decT(mas)   pmraT(mas/yr)   pmdecT(mas/yr)   e_pmraT(mas/yr)   e_pmdecT(mas/yr)  raG(deg)   e_raG(mas)   decG(deg)   e_decG(mas)   pmraG(mas/yr)   pmdecG(mas/yr)   e_pmraG(mas/yr)   e_pmdecG(mas/yr)   Parallax(mas)   e_plx(mas)   VTmag  BTmag  Teff(K)  SpClass  LClass'''

fin = open(fty2, 'r')
lin = fin.readline().strip('\n')
while len(lin):
    dat = lin.split()
    tag1= dat[0]
    ind = np.where(taggt==tag1)[0]
    if len(ind):
        j = ind[0]
        m = np.where(tagtt==tag1)[0][0]
#        print mß
        print>>fout, ('%14s   '+'%14.10f   %7.4f   '*2+'%9.6f   %9.6f   '*4)%\
        (tag1, ratt[m], ra_errtt[m], dectt[m], dec_errtt[m], \
        pmratt[m], pmra_errtt[m], pmdectt[m], pmdec_errtt[m], plxt[j], plx_errt[j],lt[j], bt[j])+\
        lin[14:]
        num += 1
    
    ind = np.where(TYC==tag1)[0]
    if len(ind):
        j = ind[0]
        hip = Hip[j]
        k = np.where(taggh==hip)[0]
        m = np.where(tagth==hip)[0][0]
#        print mß
        if len(k):
            print>>fouh, ('%14s   '+'%14.10f   %7.4f   '*2+'%10.6f   %10.6f   '*4)%\
            (hip, rath[m], ra_errth[m], decth[m], dec_errth[m], \
            pmrath[m], pmra_errth[m], pmdecth[m], pmdec_errth[m], plxh[k], plx_errh[k],lh[k], bh[k])+\
            lin[14:]
            num += 1
        
    
    N += 1
    lin = fin.readline().strip('\n')
        
print '%d '%num, 'common entries are found.'
    
print>>flog, 'Result: %d/%d entries are found. Completed!'%(num, N)
print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
flog.close()
fin.close()
fout.close()
fouh.close()

print 'Done!'