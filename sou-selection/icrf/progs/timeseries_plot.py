# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:43:55 2015

Plot coodinate time series for radio sources. 
@author: Neo
"""
import matplotlib.pyplot as plt
import fun

cat_dir = '../catalog/'
dat_dir = '../data/opa/Raw/'
res_dir = '../plot/Time_Series'

sou_list = 'sou.opa'

# read catalog file to get name of sources.
fcat = open(cat_dir + sou_list,'r')
catcon = fcat.readlines()
fcat.close()

for i in range(len(catcon)):
    sou_fil = str(catcon[i]).strip('\n')
#    print sou_fil
    fdat = open(dat_dir + sou_fil + '.txt','r')
    filcon = fdat.readlines()
    datcon = filcon[3:]
    
    t=[]
    RA=[]
    DEC=[]
    ERA=[]
    EDEC=[]

    #data format
    #Obs_moment   R.A.(o)   Dec.(o)  Err. (mas)  Corr. Del.  IERS   IVS    
    for j in range(len(datcon)-3):
        [t[j], RA[j], DEC[j], ERA[j], EDEC[j]] = \
        str(datcon[j+3]).split()[:5]
        t[j] = fun.ADepo(float(t[j]))
        RA[j] = float(RA[j]) * 3.6e9
        DEC[j] = float(DEC[j]) * 3.6e9
        ERA[j] = float(ERA[j]) * 3.6e9
        EDEC[j] = float(EDEC[j]) * 3.6e9
        
#   unit transfermation
    

    
#    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
#   ax0.errorbar(x, y, yerr=error, fmt='-o')
#    ax0.set_title('variable, symmetric error')
#
#ax1.errorbar(x, y, xerr=asymmetric_error, fmt='o')
#ax1.set_title('variable, asymmetric error')
#ax1.set_yscale('log')
#plt.show()

