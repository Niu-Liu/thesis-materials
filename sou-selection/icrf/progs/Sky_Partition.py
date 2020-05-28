# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 08:40:57 2016

@author: Neo
"""

dat_dir =  '../results/'
dat_fil = '499_sou.apm'

l1 = []
l2 = []
l3 = []
l4 = []

N1 = 30
N2 = -30
#Unit: degree
#Data are binned by intervals of declination.
#2 Nodes, +/-pi/6, dividing sources into 4 intervals.
#Each interval has the same spherical sky area.    

fdat = open(dat_dir + dat_fil, 'r')
data = fdat.readlines()

for i in range(len(data)):
    datl = data[i].strip('\n').split()
#data format
#source   P.M.(uas/yr)  Err. (uas/yr)   R.A.(o)   Dec.(o)  Err. (mas)
#            (Total R.A. Dec.)
#R.A. and Dec. corrspond to JD2000.0
    soun = datl[0]
#   mu_tol,mu_ra, mu_de, emu_tol, emu_ra, emu_de,ra0, de0, era0, ede0
    de0 = float(datl[8])
    if de0 > N1:
        l1.append(soun)
    elif de0 >0:
        l2.append(soun)
    elif de0 >N2:
        l3.append(soun)
    else:
        l4.append(soun)
        
print len(l1), len(l2), len(l3), len(l4)
