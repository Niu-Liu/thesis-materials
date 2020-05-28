# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 22:41:06 2016

@author: neo
"""

import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Apm_Load 

##for icrf1
[Sou1, V1, V_er1r, V1_RA, V1_Dec, EV_RA1, EV_Dec1, RA1, Dec1] \
= Apm_Load('icrf1_sou.apm')
##for icrf2
[Sou2, V2, V_err2, V2_RA, V2_Dec, EV_RA2, EV_Dec2, RA2, Dec2] \
= Apm_Load('icrf2_sou.apm')
##for MFV
[Sou3, V3, V_err3, V3_RA, V3_Dec, EV_RA3, EV_Dec3, RA3, Dec3] \
= Apm_Load('MFV247_sou.apm')

#set the size of subplots
left,width = 0.10,0.85
bottom_1,height = 0.1, 0.26
bottom_2 = bottom_1 + height + 0.01
bottom_3 = bottom_2 + height + 0.01
bottom_4 = bottom_3 + height + 0.01

scale3 = [left, bottom_3, width, height]
scale2 = [left, bottom_2, width, height]
scale1 = [left, bottom_1, width, height]
scale4 = [left, bottom_4, width*0.25, height/6.0]

plt.figure(figsize=(15,10))

ax1 = plt.axes(scale1)
ax2 = plt.axes(scale2)
ax3 = plt.axes(scale3, sharex = ax2, sharey = ax2)
ax4 = plt.axes(scale4, frameon=False)

m=0
for i in range(len(Sou1)):
    if V1[i]< 500:
        ax1.arrow(RA1[i], Dec1[i], V1_RA[i]/10, V1_Dec[i]/10, \
        width = 0.1, color = 'k')
    else:
        ax1.arrow(RA1[i], Dec1[i], V1_RA[i]/V1[i]*500/10, \
        V1_Dec[i]/V1[i]*500/10, \
        width = 0.1, color = 'r')
    ax1.plot(RA1[i], Dec1[i], 'b.')

    
for i in range(len(Sou2)):
    if V2[i]< 500:
        ax2.arrow(RA2[i], Dec2[i], V2_RA[i]/10, V2_Dec[i]/10, \
        width = 0.1, color = 'k')
    else:
        ax2.arrow(RA2[i], Dec2[i], V2_RA[i]/V2[i]*500/10, \
        V2_Dec[i]/V2[i]*500/10, \
        width = 0.1, color = 'r')
    ax2.plot(RA2[i], Dec2[i], 'b.')
    
for i in range(len(Sou3)):
    if V3[i]< 500:
        ax3.arrow(RA3[i], Dec3[i], V3_RA[i]/10, V3_Dec[i]/10, \
        width = 0.1, color = 'k')
    else:
        ax3.arrow(RA3[i], Dec3[i], V3_RA[i]/V3[i]*500/10, \
        V3_Dec[i]/V3[i]*500/10, \
        width = 0.1, color = 'r')
    ax3.plot(RA3[i], Dec3[i], 'b.')

ax1.set_xlim(0, 360)
ax1.set_ylim(-90, 90)
ax1.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
ax1.set_xticklabels(['0', '', '6h', '', '12h', '', '18h', '', '24h'], \
fontsize = 20)
ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax1.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 20)
#ax1.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax1.grid(True)
ax1.set_xlabel('Right Ascension', fontsize = 25)

ax2.set_xlim(0, 360)
ax2.set_ylim(-90, 90)
ax2.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
ax2.set_xticklabels(['', '', '', '', '', '', '', '', ''])
ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax2.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 20)
ax2.set_ylabel('  Declination  ($\degree$)', fontsize = 25)
ax2.grid(True)

#ax3.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax3.grid(True)
ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax3.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 20)

ax4.set_yticks([])
ax4.set_xticks([])
ax4.set_xlim(0, 90)
ax4.set_ylim(-15, 15)
ax4.plot([10, 10], [-5, 5], 'k')
ax4.plot([20, 20], [-5, 5], 'k')
ax4.plot([10, 20], [0  , 0  ], 'k')
ax4.text(25, -5, '$100\mu as/yr$', fontsize = 20)
