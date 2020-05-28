# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 22:13:30 2016

@author: Neo

Oct 25 2016: updated by Niu.
"""

import numpy as np
import matplotlib.pyplot as plt

res_dir = '../results/'
rot_fil = 'OVrot.dat'
rag_fil = 'OVrag.dat'
Num = np.loadtxt(res_dir + rot_fil, usecols=(0,), dtype=int)
W, WX, WY, WZ = \
        np.loadtxt(res_dir + rot_fil, usecols=list(range(1,9,2)), unpack=True)

w1, wx1, wy1, wz1, g1, gx1, gy1, gz1 = \
        np.loadtxt(res_dir + rag_fil, usecols=list(range(1,17,2)), unpack=True)

i = Num
z = Num
#[ z, w1, wx1, wy1, wz1, g1, gx1, gy1, gz1] = RG_Load('Rank_RG.dat')

#t1 = (np.array(W) + np.array(G))/2
#t2 = (np.array(w) + np.array(g))/2

#y = np.ones(len(i)) 
#    
#fig, ax = plt.subplots(2, 2)
#((ax1, ax2), (ax3, ax4)) = ax
##################################################
#Plot for axis rotation
#ax1.plot(i, WX, 'b')
#ax1.plot(j, wx, 'r')
#ax1.set_ylabel('$r_1$',fontsize = 25)
#ax1.plot(i, w3[0][1]*y, ':', label = '212 ICRF' )
#ax1.plot(i, w3[1][1]*y, '-.', label = '295 ICRF' )
#ax1.plot(i, w3[2][1]*y, '--', label = '247 MFV' )
#
#ax2.plot(i, WY, 'b')
#ax2.plot(j, wy, 'r')
#ax2.set_ylabel('$r_2$',fontsize = 25)
#ax2.plot(i, w3[0][2]*y, ':', label = '212 ICRF'  )
#ax2.plot(i, w3[1][2]*y, '-.', label = '295 ICRF' )
#ax2.plot(i, w3[2][2]*y, '--', label = '247 MFV' )
#
#ax3.plot(i, WZ, 'b')
#ax3.plot(j, wz, 'r')
#ax3.set_ylabel('$r_3$',fontsize = 25)
#ax3.plot(i, w3[0][3]*y, ':', label = '212 ICRF'  )
#ax3.plot(i, w3[1][3]*y, '-.', label = '295 ICRF' )
#ax3.plot(i, w3[2][3]*y, '--', label = '247 MFV' )
#    
#ax4.plot(i, W, 'b')
#ax4.plot(j, w, 'r')
#ax4.set_ylabel('$r$',fontsize = 25)
#ax4.plot(i, w3[0][0]*y, ':' , label = '212 ICRF' )
#ax4.plot(i, w3[1][0]*y, '-.', label = '295 ICRF' )
#ax4.plot(i, w3[2][0]*y, '--', label = '247 MFV' )

#####################################################
##Plot for glide
#ax1.plot(i, GX, 'b')
#ax1.plot(j, gx, 'r')
#ax1.set_ylabel('$g_1$',fontsize = 25)
#ax1.plot(i, g3[0][1]*y, ':', label = '212 ICRF' )
#ax1.plot(i, g3[1][1]*y, '-.', label = '295 ICRF' )
#ax1.plot(i, g3[2][1]*y, '--', label = '247 MFV' )
#
#ax2.plot(i, GY, 'b')
#ax2.plot(j, gy, 'r')
#ax2.set_ylabel('$g_2$',fontsize = 25)
#ax2.plot(i, g3[0][2]*y, ':', label = '212 ICRF'  )
#ax2.plot(i, g3[1][2]*y, '-.', label = '295 ICRF' )
#ax2.plot(i, g3[2][2]*y, '--', label = '247 MFV' )
#
#ax3.plot(i, GZ, 'b')
#ax3.plot(j, gz, 'r')
#ax3.set_ylabel('$g_3$',fontsize = 25)
#ax3.plot(i, g3[0][3]*y, ':', label = '212 ICRF'  )
#ax3.plot(i, g3[1][3]*y, '-.', label = '295 ICRF' )
#ax3.plot(i, g3[2][3]*y, '--', label = '247 MFV' )
#    
#ax4.plot(i, G, 'b')
#ax4.plot(j, g, 'r')
#ax4.set_ylabel('$g$',fontsize = 25)
#ax4.plot(i, g3[0][0]*y, ':' , label = '212 ICRF' )
#ax4.plot(i, g3[1][0]*y, '-.', label = '295 ICRF' )
#ax4.plot(i, g3[2][0]*y, '--', label = '247 MFV' )
##
#ax1.legend()
#ax2.legend()
#ax3.legend()
#ax4.legend()
#############################################################
## Rotation and glide, compared with rotation solely

############################################################
#common part
plt.figure()
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

ax1.plot(i, WX, 'b', linewidth=5)
ax1.plot(z, wx1, 'r.', markersize=4)
ax1.plot(z, gx1, 'g.', markersize=4)
#ax1.set_ylabel('$X$',fontsize = 25)

ax2.plot(i, WY, 'b', linewidth=5)
ax2.plot(z, wy1, 'r.', markersize=4)
ax2.plot(z, gy1, 'g.', markersize=4)
#ax2.set_ylabel('$Y$',fontsize = 25)

ax3.plot(i, WZ, 'b', linewidth=5)
ax3.plot(z, wz1, 'r.', markersize=4)
ax3.plot(z, gz1, 'g.', markersize=4)
#ax3.set_ylabel('$Z$',fontsize = 25)
    
ax4.plot(i, W, 'b', linewidth=5)
ax4.plot(z, w1, 'r.', markersize=4, label = '$r$')
ax4.plot(z, g1, 'g.', markersize=4, label = '$d$')
#ax4.set_ylabel('$\sqrt{X^2+Y^2+Z^2}$',fontsize = 25)

ax1.set_ylabel('$r_1$',fontsize = 25)
ax1.set_xlim([100, max(i)])
ax1.set_xticks([100,150,200,250,300,350,400,450,max(i)])
ax1.set_xticklabels(['','','','','','','','',''])
ax1.set_ylim([-5, 10])
ax1.set_yticks([-5,0,5])
ax1.set_yticklabels(['-5','0','5'],fontsize = 15)
#ax1.plot(i, w3[0][1]*y, ':', label = '212 ICRF' )
#ax1.plot(i, w3[1][1]*y, '-.', label = '295 ICRF' )
#ax1.plot(i, w3[2][1]*y, '--', label = '247 MFV' )

ax2.set_ylabel('$r_2$',fontsize = 25)
ax2.set_ylim([-10, 5])
ax2.set_yticks([-10,-5,0])
ax2.set_yticklabels(['-10','-5','0'],fontsize = 15)
#ax2.plot(i, w3[0][2]*y, ':', label = '212ICRF'  )
#ax2.plot(i, w3[1][2]*y, '-.', label = '295ICRF' )
#ax2.plot(i, w3[2][2]*y, '--', label = '247MFV' )

ax3.set_ylabel('$r_3$',fontsize = 25)
ax3.set_ylim([-5,10])
ax3.set_yticks([-5,0,5])
ax3.set_yticklabels(['-5','0','5'],fontsize = 15)
#ax3.plot(i, w3[0][3]*y, ':', label = '212ICRF'  )
#ax3.plot(i, w3[1][3]*y, '-.', label = '295ICRF' )
#ax3.plot(i, w3[2][3]*y, '--', label = '247MFV' )

ax4.set_ylabel('$r$',fontsize = 25)
ax4.set_ylim([0,12])
ax4.set_yticks([0,5,10])
ax4.set_yticklabels(['0','5','10'],fontsize = 15)
#ax4.plot(i, w3[0][0]*y, ':' , label = '212ICRF' )
#ax4.plot(i, w3[1][0]*y, '-.', label = '295ICRF' )
#ax4.plot(i, w3[2][0]*y, '--', label = '247MFV' )
ax4.set_xlim([100, max(i)])
ax4.set_xticks([100,150,200,250,300,350,400,450,500,550])
ax4.set_xticklabels(['100','','200','','300','','400','','500',''],fontsize = 15)
#ax4.legend(loc=2,fontsize = 10)
ax4.set_xlabel('No. Sources',fontsize = 15)
#ax4.plot(376,Wg[276]-0.5,'^')
ax4.legend(loc = 0,fontsize = 15)
plt.show()
plt.savefig('../plot/RotAndGli.eps', dpi=100)
plt.close()