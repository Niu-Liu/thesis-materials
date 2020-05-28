# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:39:38 2016

According the rank list based on linear drift, calculate the rotation of three axis

@author: Neo
"""

dat_fil = 'NC3.rank'

import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Rank_Load 
from fun import axis_rot

res_dir = '/home/neo/Documents/icrf/results/'
res = 'Rot_Rank.dat'
#res = 'Rot_Rankg.dat'
fres = open(res_dir + res,'w')

[Sou, V, V_err, V_RA,V_Dec, Va_E, Vd_E, RA, Dec] = Rank_Load(dat_fil)

W = []
W_E = []
WX = []
WX_E = []
WY = []
WY_E = []
WZ = []
WZ_E = []

w3 = [\
#w      wx       wy     wz 
#icrf1    
[17.81,    10.67,    -13.68,    4.05],\
#icrf2    
[11.31,    5.08,    -9.51,    3.43],\
#MFV247    
[10.15,    4.41,    -8.03,    4.37],\
]


for i in range(100,len(Sou)+1):
    ua = V_RA[:i+1]
    sigua = Va_E[:i+1]
    ud = V_Dec[:i+1]
    sigud = Vd_E[:i+1]
    alp0 = RA[:i+1]
    det0 = Dec[:i+1]
    
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua,ud,sigua,sigud,alp0,det0)
    
    print("%3d    %2.4f    %2.2f    %2.4f    %2.4f    %2.4f    %2.2f    %2.2f    %2.2f"%\
    (i,w,sigw,wx,wy,wz,sigwx,sigwy,sigwz), file=fres)
    
    W.append(w)
    W_E.append(sigw)
    WX.append(wx)
    WX_E.append(sigwx)
    WY.append(wy)
    WY_E.append(sigwy)
    WZ.append(wz)
    WZ_E.append(sigwz)
    
fres.close()
    
##group rank    
dat_fil = 'NC3g.rank'    
res = 'Rot_Rankg.dat'
fres = open(res_dir + res,'w')

[Sou, V, V_err, V_RA,V_Dec, Va_E, Vd_E, RA, Dec] = Rank_Load(dat_fil)

Wg = []
W_Eg = []
WXg = []
WX_Eg = []
WYg = []
WY_Eg = []
WZg = []
WZ_Eg = []

for i in range(100,len(Sou)+1):
    ua = V_RA[:i+1]
    sigua = Va_E[:i+1]
    ud = V_Dec[:i+1]
    sigud = Vd_E[:i+1]
    alp0 = RA[:i+1]
    det0 = Dec[:i+1]
    
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua,ud,sigua,sigud,alp0,det0)
    
    print("%3d    %2.4f    %2.2f    %2.4f    %2.4f    %2.4f    %2.2f    %2.2f    %2.2f"%\
    (i,w,sigw,wx,wy,wz,sigwx,sigwy,sigwz), file=fres)
    
    Wg.append(w)
    W_Eg.append(sigw)
    WXg.append(wx)
    WX_Eg.append(sigwx)
    WYg.append(wy)
    WY_Eg.append(sigwy)
    WZg.append(wz)
    WZ_Eg.append(sigwz)
    
fres.close()

#############################################################################
#plot
    
i = list(range(100,len(Sou)+1))
y = np.ones(len(i))
    
fig, ax = plt.subplots(2, 2)
((ax1, ax2), (ax3, ax4)) = ax

ax1.plot(i, WX, 'b')
ax1.plot(i, WXg, 'r')
ax1.set_ylabel('$\omega_x$',fontsize = 25)
ax1.plot(i, w3[0][1]*y, ':', label = '212 ICRF' )
ax1.plot(i, w3[1][1]*y, '-.', label = '295 ICRF' )
ax1.plot(i, w3[2][1]*y, '--', label = '247 MFV' )

ax2.plot(i, WY, 'b')
ax2.plot(i, WYg, 'r')
ax2.set_ylabel('$\omega_y$',fontsize = 25)
ax2.plot(i, w3[0][2]*y, ':', label = '212 ICRF'  )
ax2.plot(i, w3[1][2]*y, '-.', label = '295 ICRF' )
ax2.plot(i, w3[2][2]*y, '--', label = '247 MFV' )

ax3.plot(i, WZ, 'b')
ax3.plot(i, WZg, 'r')
ax3.set_ylabel('$\omega_z$',fontsize = 25)
ax3.plot(i, w3[0][3]*y, ':', label = '212 ICRF'  )
ax3.plot(i, w3[1][3]*y, '-.', label = '295 ICRF' )
ax3.plot(i, w3[2][3]*y, '--', label = '247 MFV' )
    
ax4.plot(i, W, 'b')
ax4.plot(i, Wg, 'r')
ax4.set_ylabel('$\omega$',fontsize = 25)
ax4.plot(i, w3[0][0]*y, ':' , label = '212 ICRF' )
ax4.plot(i, w3[1][0]*y, '-.', label = '295 ICRF' )
ax4.plot(i, w3[2][0]*y, '--', label = '247 MFV' )

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

ax1.set_xlabel('No. Sources',fontsize = 15)
ax2.set_xlabel('No. Sources',fontsize = 15)
ax3.set_xlabel('No. Sources',fontsize = 15)
ax4.set_xlabel('No. Sources',fontsize = 15)

plt.show()
#plt.savefig('/home/neo/Documents/icrf/plot/rot_num.eps')
    
    
print('Done!')