# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 16:17:59 2016

Consider Rotation and Glide vector both

@author: Neo
"""

#
#def Rotation_Glide(ua,ud,sigua,sigud,alp0,det0):
#        return [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,\
#            g,sigg,gx,gy,gz,siggx,siggy,siggz,chi]

from fun import Rotation_Glide
from Data_Load import Rank_Load 

res_dir = '../results/'
#dat_fil = 'NC3.rank'
#res = 'Rank_RG.dat'

#dat_fil = 'GR.rank'
#res = 'GR2.dat'
dat_fil = 'OV.rank'
res = 'OV2.dat'
fres = open(res_dir + res,'w')

[Sou, V, V_err, V_RA,V_Dec, Va_E, Vd_E, RA, Dec] = Rank_Load(dat_fil)

for i in range(100,len(Sou)+1):
    ua = V_RA[:i+1]
    sigua = Va_E[:i+1]
    ud = V_Dec[:i+1]
    sigud = Vd_E[:i+1]
    alp0 = RA[:i+1]
    det0 = Dec[:i+1]
    
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,\
            g,sigg,gx,gy,gz,siggx,siggy,siggz,chi] \
    =Rotation_Glide(ua,ud,sigua,sigud,alp0,det0)
    
    print >> fres,\
"%3d    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f    %2.4f"%\
    (i,w,wx,wy,wz,g,gx,gy,gz)
fres.close()
