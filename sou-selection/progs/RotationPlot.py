#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:38:00 2016

@author: Neo

Oct 25, 2016: updated by Niu.
"""

import numpy as np
import matplotlib.pyplot as plt

res_dir = '../results/'
dat_fil = ['OVrot.dat', 'GRrot.dat']

#glide for the three special catalog
w3 = np.loadtxt(res_dir + 'SpecialSets.rot', usecols=(1,3,5,7))

Num = np.loadtxt(res_dir + dat_fil[0], usecols=(0,), dtype=int)
W, W_E, WX, WX_E, WY, WY_E, WZ, WZ_E = \
        np.loadtxt(res_dir + dat_fil[0], usecols=list(range(1,9)), unpack=True)
Wg, W_Eg, WXg, WX_Eg, WYg, WY_Eg, WZg, WZ_Eg = \
        np.loadtxt(res_dir + dat_fil[1], usecols=list(range(1,9)), unpack=True)

i = Num
y = np.ones(Num.size)

#############################################################################
#plot, writted a year ago.
    
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

##############################################################################
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

ax1.plot(i, np.abs(WX) , 'b', linewidth=3)
ax1.plot(i, np.abs(WXg), 'r', linewidth=3)
ax1.set_ylabel('$r_1$',fontsize = 25)
ax1.set_xlim([100, max(i)])
ax1.set_xticks([100,150,200,250,300,350,400,450,500,550])
ax1.set_xticklabels(['','','','','','','','','',''])
ax1.set_ylim([0,15])
ax1.set_yticks(np.arange(0, 15, 5))
ax1.set_yticklabels(['0','5','10'],fontsize = 12)
ax1.plot(i, np.abs(w3[0][1])*y, 'b--', label = '212 ICRF' )
ax1.plot(i, np.abs(w3[1][1])*y, 'g--', label = '295 ICRF')
ax1.plot(i, np.abs(w3[2][1])*y, 'y--', label = '247 MFV' )
ax1.plot(i, np.abs(w3[3][1])*y, 'k--', label = '260 AMS' )

ax2.plot(i, np.abs(WY) , 'b', linewidth=3)
ax2.plot(i, np.abs(WYg), 'r', linewidth=3)
ax2.set_ylabel('$r_2$',fontsize = 25)
ax2.set_ylim([0,20])
ax2.set_yticks([0,5,10,15])
ax2.set_yticklabels(['0','5','10','15'],fontsize = 12)
ax2.plot(i, np.abs(w3[0][2])*y, 'b--', label = '212 ICRF' )
ax2.plot(i, np.abs(w3[1][2])*y, 'g--', label = '295 ICRF')
ax2.plot(i, np.abs(w3[2][2])*y, 'y--', label = '247 MFV' )
ax2.plot(i, np.abs(w3[3][2])*y, 'k--', label = '260 AMS' )

ax3.plot(i, np.abs(WZ) , 'b', linewidth=3)
ax3.plot(i, np.abs(WZg), 'r', linewidth=3)
ax3.set_ylabel('$r_3$',fontsize = 25)
ax3.set_ylim([0, 15])
ax3.set_yticks(np.arange(0, 15, 5))
ax3.set_yticklabels(['0','5','10'],fontsize = 12)
ax3.plot(i, w3[0][3]*y, 'b--', label = '212 ICRF' )
ax3.plot(i, w3[1][3]*y, 'g--', label = '295 ICRF')
ax3.plot(i, w3[2][3]*y, 'y--', label = '247 MFV' )
ax3.plot(i, w3[3][3]*y, 'k--', label = '260 AMS' )

ax4.plot(i, W, 'b', linewidth=3)
ax4.plot(i, Wg, 'r', linewidth=3)
ax4.set_ylabel('$r$', fontsize=25)
ax4.set_ylim([0,20])
ax4.set_yticks(np.arange(0, 20, 5))
ax4.set_yticklabels(['0','5','10','15'],fontsize = 12)
ax4.plot(i, w3[0][0]*y, 'b--' , label = '212 ICRF' )
ax4.plot(i, w3[1][0]*y, 'g--', label = '295 ICRF' )
ax4.plot(i, w3[2][0]*y, 'y--', label = '247 MFV' )
ax4.plot(i, w3[3][0]*y, 'k--', label = '260 AMS' )
ax4.set_xlim([100,max(i)])
ax4.set_xticks([100,150,200,250,300,350,400,450,500,550, max(i)])
ax4.set_xticklabels(['100','','200','','300','','400','','500', ''],fontsize = 15)
ax4.legend(loc=0, fontsize=10)
ax4.set_xlabel('No. Sources', fontsize=15)

plt.show()
plt.savefig('../plot/Rotation_No.eps', dpi=100)
plt.close()
    
print('Done!')