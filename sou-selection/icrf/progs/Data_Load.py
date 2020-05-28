# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 08:35:14 2016

Load the catalog and data files of ICRF

@author: Neo
"""
dat_dir =  '../results/'

def Apm_Load(dat_fil):
#    dat_dir =  '../results/'
    fdat = open(dat_dir + dat_fil, 'r')
    data = fdat.readlines()

    soun = []    
    
    mu_tol = []
    mu_ra  = []
    mu_de  = []
    
    emu_tol = []
    emu_ra  = []
    emu_de  = []
    
    ra0 = []
    de0 = []
    
#    era0 = []
#    ede0 = []

    for i in range(len(data)):
        datl = data[i].strip('\n').split()
#data format
#source   P.M.(uas/yr)  Err. (uas/yr)   R.A.(o)   Dec.(o)  Err. (mas)
#            (Total R.A. Dec.)
#R.A. and Dec. corrspond to JD2000.0
        soun.append(datl[0])
#   mu_tol,mu_ra, mu_de, emu_tol, emu_ra, emu_de,ra0, de0, era0, ede0
        mu_tol.append(float(datl[1]))
        mu_ra.append(float(datl[2]))
        mu_de.append(float(datl[3]))

        emu_tol.append(float(datl[4]))
        emu_ra.append(float(datl[5]))
        emu_de.append(float(datl[6]))
    
        ra0.append(float(datl[7]))
        de0.append(float(datl[8]))
        
#        era0.append(float(datl[9]))
#        ede0.append(float(datl[10]))
        
    return [ \
        soun,\
        mu_tol,emu_tol, \
        mu_ra, mu_de, \
        emu_ra, emu_de,\
        ra0, de0 \
#        era0, ede0 \
        ]
        
def Rank_Load(dat_fil):
#    dat_dir =  '../results/'
    fdat = open(dat_dir + dat_fil, 'r')
    data = fdat.readlines()

    soun = []    
    
    mu_tol = []
    mu_ra  = []
    mu_de  = []
    
    emu_tol = []
    emu_ra  = []
    emu_de  = []
    
    ra0 = []
    de0 = []

    for i in range(len(data)):
        datl = data[i].strip('\n').split()
#data format
#        [('sou', 'S10'), ('u', float), ('u_e', float),\
#        ('alp', float), ('ua', float), ('ua_e', float),\
#        ('det', float), ('ud', float), ('ud_e', float)]
#R.A. and Dec. corrspond to JD2000.0
        soun.append(datl[0])
#   mu_tol,mu_ra, mu_de, emu_tol, emu_ra, emu_de,ra0, de0, era0, ede0
        mu_tol.append(float(datl[1]))
        mu_ra.append(float(datl[4]))
        mu_de.append(float(datl[7]))

        emu_tol.append(float(datl[2]))
        emu_ra.append(float(datl[5]))
        emu_de.append(float(datl[8]))
    
        ra0.append(float(datl[3]))
        de0.append(float(datl[6]))
        
    return [ soun,\
        mu_tol,emu_tol,\
        mu_ra, mu_de, \
        emu_ra, emu_de,\
        ra0, de0 ]
        
def Rot_Load(dat_fil):
#    dat_dir =  '../results/'
    fdat = open(dat_dir + dat_fil, 'r')
    data = fdat.readlines()

    ind = []    
    
    w = []
    wx  = []
    wy  = []
    wz = []

    for i in range(len(data)):
        datl = data[i].strip('\n').split()
        ind.append(int(datl[0]))
        w.append(float(datl[1]))
        wx.append(float(datl[3]))
        wy.append(float(datl[4]))
        wz.append(float(datl[5]))
        
    return [ ind, w, wx, wy, wz ]
    
def RG_Load(dat_fil):
#    dat_dir =  '../results/'
    fdat = open(dat_dir + dat_fil, 'r')
    data = fdat.readlines()

    ind = []    
    
    w = []
    wx  = []
    wy  = []
    wz = []
    
    g = []
    gx  = []
    gy  = []
    gz = []

    for i in range(len(data)):
        datl = data[i].strip('\n').split()
        ind.append(int(datl[0]))
        
        w.append(float(datl[1]))
        wx.append(float(datl[2]))
        wy.append(float(datl[3]))
        wz.append(float(datl[4]))
        
        g.append(float(datl[5]))
        gx.append(float(datl[6]))
        gy.append(float(datl[7]))
        gz.append(float(datl[8]))
        
    return [ ind, w, wx, wy, wz, g, gx, gy, gz ]    
    
