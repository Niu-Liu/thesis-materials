# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 19:34:26 2016

@author: neo
"""

import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Apm_Load 

[Sou1, V1, V_er1r, V1_RA, V1_Dec, EV_RA1, EV_Dec1, RA1, Dec1] \
= Apm_Load('322_sou.apm')
[Sou2, V2, V_err2, V2_RA, V2_Dec, EV_RA2, EV_Dec2, RA2, Dec2] \
= Apm_Load('351_sou.apm')
[Sou3, V3, V_err3, V3_RA, V3_Dec, EV_RA3, EV_Dec3, RA3, Dec3] \
= Apm_Load('371_sou.apm')
[Sou4, V4, V_err4, V4_RA, V4_Dec, EV_RA4, EV_Dec4, RA4, Dec4] \
= Apm_Load('394_sou.apm')

#ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

#set the size of subplots
left,width = 0.10,0.85
bottom_1,height = 0.05, 0.20
bottom_2 = bottom_1 + height + 0.01
bottom_3 = bottom_2 + height + 0.01
bottom_4 = bottom_3 + height + 0.01
bottom_5 = bottom_4 + height + 0.01

scale3 = [left, bottom_3, width, height]
scale2 = [left, bottom_2, width, height]
scale1 = [left, bottom_1, width, height]
scale4 = [left, bottom_4, width, height]
scale5 = [left, bottom_5, width*0.25, height/6.0]

plt.figure(figsize=(50,10))

ax1 = plt.axes(scale1)
ax2 = plt.axes(scale2)
ax3 = plt.axes(scale3, sharex = ax2, sharey = ax2)
ax4 = plt.axes(scale4, sharex = ax2, sharey = ax2)
ax5 = plt.axes(scale5, frameon=False)

for i in range(len(Sou1)):
#    if V1[i]< 500:
    ax1.arrow(RA1[i], Dec1[i], V1_RA[i]/10, V1_Dec[i]/10, \
        width = 0.1, color = 'k')
#    else:
#        ax1.arrow(RA1[i], Dec1[i], V1_RA[i]/V1[i]*500/10, \
#        V1_Dec[i]/V1[i]*500/10, \
#        width = 0.1, color = 'r')
#    ax1.plot(RA1[i], Dec1[i], 'b.')
    if Sou1[i] in defn:
        ax1.plot(RA1[i], Dec1[i], 'ro')
    else:
        ax1.plot(RA1[i], Dec1[i], 'b.')

    
for i in range(len(Sou2)):
#    if V2[i]< 500:
    ax2.arrow(RA2[i], Dec2[i], V2_RA[i]/10, V2_Dec[i]/10, \
        width = 0.1, color = 'k')
#    else:
#        ax2.arrow(RA2[i], Dec2[i], V2_RA[i]/V2[i]*500/10, \
#        V2_Dec[i]/V2[i]*500/10, \
#        width = 0.1, color = 'r')
#    ax2.plot(RA2[i], Dec2[i], 'b.')
    if Sou2[i] in defn:
        ax2.plot(RA2[i], Dec2[i], 'ro')
    else:
        ax2.plot(RA2[i], Dec2[i], 'b.')    
    
for i in range(len(Sou3)):
#    if V3[i]< 500:
    ax3.arrow(RA3[i], Dec3[i], V3_RA[i]/10, V3_Dec[i]/10, \
        width = 0.1, color = 'k')
#    else:
#        ax3.arrow(RA3[i], Dec3[i], V3_RA[i]/V3[i]*500/10, \
#        V3_Dec[i]/V3[i]*500/10, \
#        width = 0.1, color = 'r')
#    ax3.plot(RA3[i], Dec3[i], 'b.')
    if Sou3[i] in defn:
        ax3.plot(RA3[i], Dec3[i], 'ro')
    else:
        ax3.plot(RA3[i], Dec3[i], 'b.')

for i in range(len(Sou4)):
    ax4.arrow(RA4[i], Dec4[i], V4_RA[i]/10, V4_Dec[i]/10, \
        width = 0.1, color = 'k')
#    ax4.plot(RA4[i], Dec4[i], 'b.')
    if Sou4[i] in defn:
        ax4.plot(RA4[i], Dec4[i], 'ro')
    else:
        ax4.plot(RA4[i], Dec4[i], 'b.')    

ax1.set_xlim(0, 360)
ax1.set_ylim(-90, 90)
ax1.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
ax1.set_xticklabels(['0', '', '6h', '', '12h', '', '18h', '', '24h'], \
fontsize = 15)
ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax1.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 15)
#ax1.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax1.grid(True)
ax1.set_xlabel('Right Ascension', fontsize = 15)
ax1.set_ylabel('  Declination  ($\degree$)', fontsize = 15)

ax2.set_xlim(0, 360)
ax2.set_ylim(-90, 90)
ax2.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
ax2.set_xticklabels(['', '', '', '', '', '', '', '', ''])
ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax2.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 15)
ax2.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax2.grid(True)

#ax3.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax3.grid(True)
ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax3.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 15)
ax3.set_ylabel('  Declination  ($\degree$)', fontsize = 15)

ax4.grid(True)
ax4.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax4.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 15)
ax4.set_ylabel('  Declination  ($\degree$)', fontsize = 15)

ax5.set_yticks([])
ax5.set_xticks([])
ax5.set_xlim(0, 90)
ax5.set_ylim(-15, 15)
ax5.plot([10, 10], [-5, 5], 'k')
ax5.plot([20, 20], [-5, 5], 'k')
ax5.plot([10, 20], [0  , 0  ], 'k')
ax5.text(25, -5, '$100\mu as/yr$', fontsize = 15)