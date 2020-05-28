# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 22:28:32 2016

For the specific sets of sources, Calculate their axis rotation.

@author: Neo
"""

dat_fil = ['icrf1', 'icrf2', 'MFV247', 'AMS260']

res_dir = '../results/'
resr = '4sets.rot'
resg = '4sets.gli'
resgr = '4sets.gr'
fresr = open(res_dir + resr,'w')
fresg = open(res_dir + resg,'w')
fresgr = open(res_dir + resgr,'w')

from Data_Load import Apm_Load 
from fun import axis_rot, simulation, Rotation_Glide
import numpy as np

for i in range(len(dat_fil)):
    dat = dat_fil[i] + '_sou.apm'
    [Sou, V, V_err, ua, ud, sigua, sigud, alp0, det0] = Apm_Load(dat)
    
##Rotation    
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua,ud,sigua,sigud,alp0,det0)
    
    print >> fresr, dat_fil[i] + \
"    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f"%\
    (w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)
    
##Rotation and Glide
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,\
            g,sigg,gx,gy,gz,siggx,siggy,siggz,chi] \
    =Rotation_Glide(ua,ud,sigua,sigud,alp0,det0)
    print >> fresgr, dat_fil[i] + \
"    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f"%\
    (w,wx,wy,wz,g,gx,gy,gz)
    
#Glide(Simulation)
    [uag, udg] = simulation(alp0, det0)
    sig = np.ones(len(uag))
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(uag, udg, sig, sig, alp0, det0)
    
    print >> fresg, dat_fil[i] + \
"    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f"%\
    (w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)
      
fresr.close()  
fresg.close() 
fresgr.close() 
print "Done!" 
