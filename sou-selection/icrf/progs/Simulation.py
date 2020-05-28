# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 16:47:01 2016

@author: Neo
"""


####################################################################
##Test the codes
#from fun import simulation, axis_rot
#from Data_Load import Apm_Load 
#import numpy as np
#dat_fil = 'icrf2_sou.apm'
#[Sou, V, V_err, V_RA, V_Dec, Va_E, Vd_E, RA, Dec] = Apm_Load(dat_fil)
#[ua, ud] = simulation(RA, Dec)
#print ua[0],ud[0]
#sig = np.ones(len(ua))
#
#[w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
#    =axis_rot(ua, ud, sig, sig, RA, Dec)
#    
#print w,wx,wy,wz
##Test finished
#####################################################################

from fun import simulation, axis_rot
from Data_Load import Rank_Load 
import numpy as np
import matplotlib.pyplot as plt

#glide for the three special catalog
w3 = [\
#icrf1    
[0.67,    0.63,    0.24,    -0.07] ,\
#icrf2    
[0.33,    0.14,    0.30,    -0.00] ,\
#MFV247    
[0.69,    -0.07,    0.69,    -0.01]  \
]

res_dir = '../results/'
dat_fil = ['OV.rank', 'GR.rank']
res = ['OVsim.dat', 'GRsim.dat']

W = []
W_E = []
WX = []
WX_E = []
WY = []
WY_E = []
WZ = []
WZ_E = []

for d in range(len(dat_fil)):
    fres = open(res_dir + res[d],'w')
    [Sou, V, V_err, V_RA,V_Dec, Va_E, Vd_E, RA, Dec] = Rank_Load(dat_fil[d])
    [ua, ud] = simulation(RA, Dec)
    sig = np.ones(len(ua))
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua, ud, sig, sig, RA, Dec)
    
    N= len(Sou)
    
    for i in range(100,len(Sou)+1):
        sig = np.ones(i)
        [ua, ud] = simulation(RA[:i], Dec[:i])
        
        [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua, ud, sig, sig, RA[:i], Dec[:i])
    
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
    
N = N-100+1    
    
Wg = W[N:]
W_Eg = W_E[N:]
WXg = WX[N:]
WX_Eg = WX_E[N:]
WYg = WY[N:]
WY_Eg = WY_E[N:]
WZg = WZ[N:]
WZ_Eg = WZ_E[N:]

W = W[:N]
W_E = W_E[:N]
WX = WX[:N]
WX_E = WX_E[:N]
WY = WY[:N]
WY_E = WY_E[:N]
WZ = WZ[:N]
WZ_E = WZ_E[:N]

i = range(100,len(Sou)+1)
y = np.ones(len(i))
    
#fig, ax = plt.subplots(2, 2)
#((ax1, ax2), (ax3, ax4)) = ax
#
#ax1.plot(i, WX, 'b')
#ax1.plot(i, WXg, 'r')
#ax1.set_ylabel('$g_x$',fontsize = 25)
#ax1.plot(i, w3[0][1]*y, ':', label = '212 ICRF' )
#ax1.plot(i, w3[1][1]*y, '-.', label = '295 ICRF' )
#ax1.plot(i, w3[2][1]*y, '--', label = '247 MFV' )
#
#ax2.plot(i, WY, 'b')
#ax2.plot(i, WYg, 'r')
#ax2.set_ylabel('$g_y$',fontsize = 25)
#ax2.plot(i, w3[0][2]*y, ':', label = '212 ICRF'  )
#ax2.plot(i, w3[1][2]*y, '-.', label = '295 ICRF' )
#ax2.plot(i, w3[2][2]*y, '--', label = '247 MFV' )
#
#ax3.plot(i, WZ, 'b')
#ax3.plot(i, WZg, 'r')
#ax3.set_ylabel('$g_z$',fontsize = 25)
#ax3.plot(i, w3[0][3]*y, ':', label = '212 ICRF'  )
#ax3.plot(i, w3[1][3]*y, '-.', label = '295 ICRF' )
#ax3.plot(i, w3[2][3]*y, '--', label = '247 MFV' )
#    
#ax4.plot(i, W, 'b')
#ax4.plot(i, Wg, 'r')
#ax4.set_ylabel('$g$',fontsize = 25)
#ax4.plot(i, w3[0][0]*y, ':' , label = '212 ICRF' )
#ax4.plot(i, w3[1][0]*y, '-.', label = '295 ICRF' )
#ax4.plot(i, w3[2][0]*y, '--', label = '247 MFV' )
#
#ax1.legend()
#ax2.legend()
#ax3.legend(loc=4)
#ax4.legend()
#
#ax1.set_xlabel('No. Sources',fontsize = 15)
#ax2.set_xlabel('No. Sources',fontsize = 15)
#ax3.set_xlabel('No. Sources',fontsize = 15)
#ax4.set_xlabel('No. Sources',fontsize = 15)
    
##Acutually only plot for g is needed.
plt.plot(i, W, 'b', linewidth=3)
plt.plot(i, Wg, 'r', linewidth=3)
plt.ylabel('$g$',fontsize = 30)
plt.plot(i, w3[0][0]*y, ':' , label = '212 ICRF' )
plt.plot(i, w3[1][0]*y, '-.', label = '295 ICRF' )
plt.plot(i, w3[2][0]*y, '--', label = '247 MFV' )
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 20)
plt.xlim([100, max(i)])
plt.legend()
plt.xlabel('No. Sources',fontsize = 15)

plt.show()
