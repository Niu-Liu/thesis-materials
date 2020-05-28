# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:04:15 2016

from Every group pick the equal number of source

@author: Neo
"""

from Data_Load import Rank_Load
import numpy as np
from fun import axis_rot, simulation
import matplotlib.pyplot as plt

dat_fil = 'NC3.rank'
res_dir = '../results/'
res = 'Eq_Rank.rot'
fres = open(res_dir + res,'w')
resg = 'Eq_Rank.gli'
fresg = open(res_dir + resg,'w')

[Sou, V, V_err, V_RA,V_Dec, Va_E, Vd_E, RA, Dec] = Rank_Load(dat_fil)

#Num = range(100/4, 496/4, 1)
dtype = [('sou', 'S10'), ('u', float), ('u_e', float),\
        ('alp', float), ('ua', float), ('ua_e', float),\
        ('det', float), ('ud', float), ('ud_e', float)]
        
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
#icrf1    \
[17.81,    10.67,    -13.68,    4.05],\
#icrf2    
[11.31,    5.08,    -9.51,    3.43],\
#MFV247    
[10.15,    4.41,    -8.03,    4.37],\
]

num = []

nodes = [-90, -30, 0, 30, 90]
for Num in range(25, 198+1):
    dat = []
    for j in range(len(nodes)-1):
        a = []
        for i in range(len(V)):
            if nodes[j] < Dec[i] < nodes[j+1]:
                a.append((Sou[i], V[i], V_err[i],\
        RA[i], V_RA[i], Va_E[i], Dec[i], V_Dec[i],Vd_E[i]))
        
        a = np.array(a, dtype=dtype)        
        a = np.sort(a, order=['u']) 
        a = list(a[: Num])
       
        dat += a

    ua = [ dat[i][4] for i in range(len(dat))]
    sigua = [ dat[i][5] for i in range(len(dat))]
    ud = [ dat[i][7] for i in range(len(dat))]
    sigud = [ dat[i][8] for i in range(len(dat))]
    alp0 = [ dat[i][3] for i in range(len(dat))]
    det0 = [ dat[i][6] for i in range(len(dat))]
    
    n = len(ua)
    
    num.append(len(ua))
    
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua,ud,sigua,sigud,alp0,det0)
    
    print >> fres,\
"%3d    %2.4f    %2.2f    %2.4f    %2.4f    %2.4f    %2.2f    %2.2f    %2.2f"%\
    (n,w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)  
    
    sig = np.ones(n)    
    [ua, ud] = simulation(alp0,det0)
    [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]\
    =axis_rot(ua,ud,sig,sig,alp0,det0)
    print >> fresg,\
"%3d    %2.4f    %2.2f    %2.4f    %2.4f    %2.4f    %2.2f    %2.2f    %2.2f"%\
    (n,w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)
    
#    W.append(w)
#    W_E.append(sigw)
#    WX.append(wx)
#    WX_E.append(sigwx)
#    WY.append(wy)
#    WY_E.append(sigwy)
#    WZ.append(wz)
#    WZ_E.append(sigwz)
#    
#i = num
#y = np.ones(len(i))
#    
#fig, ax = plt.subplots(2, 2)
#((ax1, ax2), (ax3, ax4)) = ax
#
#ax1.plot(i, WX, 'b')
#ax1.set_ylabel('$\omega_x$',fontsize = 25)
#ax1.plot(i, w3[0][1]*y, ':', label = '212 ICRF' )
#ax1.plot(i, w3[1][1]*y, '-.', label = '295 ICRF' )
#ax1.plot(i, w3[2][1]*y, '--', label = '247 MFV' )
#
#ax2.plot(i, WY, 'b')
#ax2.set_ylabel('$\omega_y$',fontsize = 25)
#ax2.plot(i, w3[0][2]*y, ':', label = '212 ICRF'  )
#ax2.plot(i, w3[1][2]*y, '-.', label = '295 ICRF' )
#ax2.plot(i, w3[2][2]*y, '--', label = '247 MFV' )
#
#ax3.plot(i, WZ, 'b')
#ax3.set_ylabel('$\omega_z$',fontsize = 25)
#ax3.plot(i, w3[0][3]*y, ':', label = '212 ICRF'  )
#ax3.plot(i, w3[1][3]*y, '-.', label = '295 ICRF' )
#ax3.plot(i, w3[2][3]*y, '--', label = '247 MFV' )
#    
#ax4.plot(i, W, 'b')
#ax4.set_ylabel('$\omega$',fontsize = 25)
#ax4.plot(i, w3[0][0]*y, ':' , label = '212 ICRF' )
#ax4.plot(i, w3[1][0]*y, '-.', label = '295 ICRF' )
#ax4.plot(i, w3[2][0]*y, '--', label = '247 MFV' )
#
#ax1.legend()
#ax2.legend(loc=4)
#ax3.legend(loc=4)
#ax4.legend()
#
#ax1.set_xlabel('No. Sources',fontsize = 15)
#ax2.set_xlabel('No. Sources',fontsize = 15)
#ax3.set_xlabel('No. Sources',fontsize = 15)
#ax4.set_xlabel('No. Sources',fontsize = 15)
#
#plt.show()
##plt.savefig('../plot/rot_num.eps')
#    
fres.close()
fresg.close()
#    
#print 'Done!'    
