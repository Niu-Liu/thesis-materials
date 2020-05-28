# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:39:38 2016

According the rank list based on linear drift, calculate the rotation of three axis

@author: Neo
"""

dat_fil = 'OV.rank'


import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Rank_Load 
from fun import axis_rot
w3 = [\
#w      wx       wy     wz 
#icrf1    
[17.81,    -10.67,    13.68,    -4.05],\
#icrf2    
[11.31,    -5.08,    9.51,    -3.43],\
#MFV247    
[10.15,    -4.41,    8.03,    -4.37],\
]

res_dir = '../results/'
res = 'OV.dat'
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

for i in range(100,len(Sou)+1):
    ua = V_RA[:i+1]
    sigua = Va_E[:i+1]
    ud = V_Dec[:i+1]
    sigud = Vd_E[:i+1]
    alp0 = RA[:i+1]
    det0 = Dec[:i+1]
    
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua,ud,sigua,sigud,alp0,det0)
    
    print >> fres,\
"%3d    %2.4f    %2.2f    %2.4f    %2.4f    %2.4f    %2.2f    %2.2f    %2.2f"%\
    (i,w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)
    
    W.append(w)
    W_E.append(sigw)
    WX.append(wx)
    WX_E.append(sigwx)
    WY.append(wy)
    WY_E.append(sigwy)
    WZ.append(wz)
    WZ_E.append(sigwz)
    
fres.close()
    
#group rank    
dat_fil = 'GR.rank'    
res = 'GR.dat'
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
    
    print >> fres,\
"%3d    %2.4f    %2.2f    %2.4f    %2.4f    %2.4f    %2.2f    %2.2f    %2.2f"%\
    (i,w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)
    
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
    
#i = range(100,len(Sou)+1)
#y = np.ones(len(i))
#    
#fig, ax = plt.subplots(2, 2)
#((ax1, ax2), (ax3, ax4)) = ax
#
#ax1.plot(i, WX, 'b')
##ax1.plot(i, WXg, 'r')
#ax1.set_ylabel('$\omega_x$',fontsize = 25)
#ax1.plot(i, w3[0][1]*y, ':', label = '212 ICRF' )
#ax1.plot(i, w3[1][1]*y, '-.', label = '295 ICRF' )
#ax1.plot(i, w3[2][1]*y, '--', label = '247 MFV' )
#
#ax2.plot(i, WY, 'b')
##ax2.plot(i, WYg, 'r')
#ax2.set_ylabel('$\omega_y$',fontsize = 25)
#ax2.plot(i, w3[0][2]*y, ':', label = '212 ICRF'  )
#ax2.plot(i, w3[1][2]*y, '-.', label = '295 ICRF' )
#ax2.plot(i, w3[2][2]*y, '--', label = '247 MFV' )
#
#ax3.plot(i, WZ, 'b')
##ax3.plot(i, WZg, 'r')
#ax3.set_ylabel('$\omega_z$',fontsize = 25)
#ax3.plot(i, w3[0][3]*y, ':', label = '212 ICRF'  )
#ax3.plot(i, w3[1][3]*y, '-.', label = '295 ICRF' )
#ax3.plot(i, w3[2][3]*y, '--', label = '247 MFV' )
#    
#ax4.plot(i, W, 'b')
##ax4.plot(i, Wg, 'r')
#ax4.set_ylabel('$\omega$',fontsize = 25)
#ax4.plot(i, w3[0][0]*y, ':' , label = '212 ICRF' )
#ax4.plot(i, w3[1][0]*y, '-.', label = '295 ICRF' )
#ax4.plot(i, w3[2][0]*y, '--', label = '247 MFV' )
#
#ax1.legend()
#ax2.legend()
#ax3.legend()
#ax4.legend()
#
#ax1.set_xlabel('No. Sources',fontsize = 15)
#ax2.set_xlabel('No. Sources',fontsize = 15)
#ax3.set_xlabel('No. Sources',fontsize = 15)
#ax4.set_xlabel('No. Sources',fontsize = 15)
#
#plt.show()
##plt.savefig('../plot/OARank_rot.eps')
#plt.savefig('../plot/GRank_rot.eps') 

## Another plot   
i = range(100,len(Sou)+1)
y = np.ones(len(i))
plt.figure()

#set the size of subplots
left,width = 0.10,0.85
bottom,height = 0.1, 0.17
bottom_3 = bottom + height*2 + 0.01
bottom_2 = bottom_3 + height + 0.01
bottom_1 = bottom_2 + height + 0.01

scale4 = [left, bottom,   width, height*2]
scale3 = [left, bottom_3, width, height]
scale2 = [left, bottom_2, width, height]
scale1 = [left, bottom_1, width, height]

ax1 = plt.axes(scale1)
ax2 = plt.axes(scale2, sharex = ax1)
ax3 = plt.axes(scale3, sharex = ax1)
ax4 = plt.axes(scale4)

ax1.plot(i, WX, 'b', linewidth=3)
ax1.plot(i, WXg, 'r', linewidth=3)
ax1.set_ylabel('$r_1$',fontsize = 25)
ax1.set_xlim([100,max(i)])
ax1.set_xticks([100,150,200,250,300,350,400,450,500,550])
ax1.set_xticklabels(['','','','','','','','','',''])
ax1.set_ylim([-15,5])
ax1.set_yticks([-15,-10,-5,0,5])
ax1.set_yticklabels(['','-10','-5','0',''],fontsize = 15)
ax1.plot(i, w3[0][1]*y, ':', label = '212 ICRF' )
ax1.plot(i, w3[1][1]*y, '-.', label = '295 ICRF' )
ax1.plot(i, w3[2][1]*y, '--', label = '247 MFV' )

ax2.plot(i, WY, 'b', linewidth=3)
ax2.plot(i, WYg, 'r', linewidth=3)
ax2.set_ylabel('$r_2$',fontsize = 25)
ax2.set_ylim([-5,15])
ax2.set_yticks([-5,0,5,10,15])
ax2.set_yticklabels(['','0','5','10',''],fontsize = 15)
ax2.plot(i, w3[0][2]*y, ':', label = '212ICRF'  )
ax2.plot(i, w3[1][2]*y, '-.', label = '295ICRF' )
ax2.plot(i, w3[2][2]*y, '--', label = '247MFV' )

ax3.plot(i, WZ, 'b', linewidth=3)
ax3.plot(i, WZg, 'r', linewidth=3)
ax3.set_ylabel('$r_3$',fontsize = 25)
ax3.set_ylim([-5,-1])
ax3.set_yticks([-5,-4,-3,-2,-1])
ax3.set_yticklabels(['','-4','-3','-2',''],fontsize = 15)
ax3.plot(i, w3[0][3]*y, ':', label = '212ICRF'  )
ax3.plot(i, w3[1][3]*y, '-.', label = '295ICRF' )
ax3.plot(i, w3[2][3]*y, '--', label = '247MFV' )

ax4.plot(i, W, 'b', linewidth=3)
ax4.plot(i, Wg, 'r', linewidth=3)
ax4.set_ylabel('$r$',fontsize = 25)
ax4.set_ylim([0,20])
ax4.set_yticks([0,5,10,15,20])
ax4.set_yticklabels(['0','5','10','15',''],fontsize = 15)
ax4.plot(i, w3[0][0]*y, ':' , label = '212ICRF' )
ax4.plot(i, w3[1][0]*y, '-.', label = '295ICRF' )
ax4.plot(i, w3[2][0]*y, '--', label = '247MFV' )
ax4.set_xlim([100,max(i)])
ax4.set_xticks([100,150,200,250,300,350,400,450,500,550])
ax4.set_xticklabels(['100','','200','','300','','400','','500',''],fontsize = 15)
ax4.legend(loc=2,fontsize = 15)
ax4.set_xlabel('No. Sources',fontsize = 15)

plt.show()
    
print 'Done!'