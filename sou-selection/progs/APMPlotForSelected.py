# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 19:34:26 2016

@author: neo

Oct 25, 2016: updated by Niu.
"""

import numpy as np
#import matplotlib as mpl
#mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

res_dir = '../results/'
dat_fil = ['OV.rank', 'GR.rank']

Souo = np.loadtxt(res_dir + dat_fil[0], usecols=(0,), dtype=str)
RAo, V_RAo, Deco, V_Deco = \
        np.loadtxt(res_dir + dat_fil[0], usecols=[3,4,6,7], unpack=True)
Soug = np.loadtxt(res_dir + dat_fil[1], usecols=(0,), dtype=str)
RAg, V_RAg, Decg, V_Decg = \
        np.loadtxt(res_dir + dat_fil[1], usecols=[3,4,6,7], unpack=True)
        
## For Our selected subsets.
ind = np.array([120, 223, 130, 194])+100
#ind = np.array([223, 368, 194, 347])+100
Sou1, V1_RA, V1_Dec, RA1, Dec1 \
= Souo[:ind[0]+1], V_Deco[:ind[0]+1], V_RAo[:ind[0]+1], RAo[:ind[0]+1], Deco[:ind[0]+1] 
Sou2, V2_RA, V2_Dec, RA2, Dec2 \
= Souo[:ind[1]+1], V_Deco[:ind[1]+1], V_RAo[:ind[1]+1], RAo[:ind[1]+1], Deco[:ind[1]+1] 
Sou3, V3_RA, V3_Dec, RA3, Dec3 \
= Soug[:ind[2]+1], V_Decg[:ind[2]+1], V_RAg[:ind[2]+1], RAg[:ind[2]+1], Decg[:ind[2]+1] 
Sou4, V4_RA, V4_Dec, RA4, Dec4 \
= Soug[:ind[3]+1], V_Decg[:ind[3]+1], V_RAg[:ind[3]+1], RAg[:ind[3]+1], Decg[:ind[3]+1]

#ICRF2 defining sources data file
icrf2f = '../data/icrf2-defining.dat'
defn = np.loadtxt(icrf2f, dtype=str, skiprows=20, usecols=(2,))

#set the size of subplots
left,width = 0.10,0.85
bottom_1,height = 0.06, 0.21
bottom_2 = bottom_1 + height + 0.01
bottom_3 = bottom_2 + height + 0.01
bottom_4 = bottom_3 + height + 0.01
bottom_5 = bottom_4 + height + 0.01

scale1 = [left, bottom_1, width, height]
scale2 = [left, bottom_2, width, height]
scale3 = [left, bottom_3, width, height]
scale4 = [left, bottom_4, width, height]
scale5 = [left, bottom_5, width*0.25, height/6.0]

plt.figure(figsize=(9, 28))

ax1 = plt.axes(scale1)
ax2 = plt.axes(scale2)
ax3 = plt.axes(scale3, sharex = ax2, sharey = ax2)
ax4 = plt.axes(scale4, sharex = ax2, sharey = ax2)
ax5 = plt.axes(scale5, frameon=False)

for i in range(len(Sou1)):
    if Sou1[i] in defn:
        ax1.plot(RA1[i], Dec1[i], 'bo', markersize=5)
        ax1.arrow(RA1[i], Dec1[i], V1_RA[i]/10, V1_Dec[i]/10, \
                  width = 0.15, color = 'b')
    else:
        ax1.plot(RA1[i], Dec1[i], 'go', markersize=5)
        ax1.arrow(RA1[i], Dec1[i], V1_RA[i]/10, V1_Dec[i]/10, \
        width = 0.15, color = 'g')
    
for i in range(len(Sou2)):
    if Sou2[i] in defn:
        ax2.plot(RA2[i], Dec2[i], 'bo', markersize=5)
        ax2.arrow(RA2[i], Dec2[i], V2_RA[i]/10, V2_Dec[i]/10, \
                  width = 0.15, color = 'b')
    else:
        ax2.plot(RA2[i], Dec2[i], 'go', markersize=5)    
        ax2.arrow(RA2[i], Dec2[i], V2_RA[i]/10, V2_Dec[i]/10, \
                  width = 0.15, color = 'g')
        
for i in range(len(Sou3)):
    if Sou3[i] in defn:
        ax3.plot(RA3[i], Dec3[i], 'bo', markersize=5)
        ax3.arrow(RA3[i], Dec3[i], V3_RA[i]/10, V3_Dec[i]/10, \
                  width = 0.15, color = 'b')
    else:
        ax3.plot(RA3[i], Dec3[i], 'go', markersize=5)
        ax3.arrow(RA3[i], Dec3[i], V3_RA[i]/10, V3_Dec[i]/10, \
                  width = 0.15, color = 'g')

for i in range(len(Sou4)):
    if Sou4[i] in defn:
        ax4.plot(RA4[i], Dec4[i], 'bo', markersize=5)
        ax4.arrow(RA4[i], Dec4[i], V4_RA[i]/10, V4_Dec[i]/10, \
                  width = 0.15, color = 'b')
    else:
        ax4.plot(RA4[i], Dec4[i], 'go', markersize=5)    
        ax4.arrow(RA4[i], Dec4[i], V4_RA[i]/10, V4_Dec[i]/10, \
                  width = 0.15, color = 'g')
    
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
        
ax1.set_xlim(0, 360)
ax1.set_ylim(-90, 90)
ax1.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
ax1.set_xticklabels(['0', '', '6h', '', '12h', '', '18h', '', '24h'], \
fontsize = 12)
ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax1.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 12)
ax1.grid(True)
ax1.set_xlabel('Right Ascension', fontsize = 15)
ax1.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax1.text( 5, 72, r'Sou%d'%ind[0], fontdict=font)

ax2.set_xlim(0, 360)
ax2.set_ylim(-90, 90)
ax2.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
ax2.set_xticklabels(['', '', '', '', '', '', '', '', ''])
ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax2.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 12)
ax2.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax2.grid(True)
ax2.text( 5, 72, r'Sou%d'%ind[1], fontdict=font)

ax3.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax3.grid(True)
ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax3.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 12)
ax3.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax3.text( 5, 72, r'Sou%d'%ind[2], fontdict=font)

ax4.grid(True)
ax4.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax4.set_yticklabels(['', '-60', '-30', '0', '30', '60', ''], fontsize = 12)
ax4.set_ylabel('  Declination  ($\degree$)', fontsize = 15)
ax4.text( 5, 72, r'Sou%d'%ind[3], fontdict=font)

ax5.set_yticks([])
ax5.set_xticks([])
ax5.set_xlim(0, 90)
ax5.set_ylim(-15, 15)
ax5.plot([10, 10], [-5, 5], 'k')
ax5.plot([20, 20], [-5, 5], 'k')
ax5.plot([10, 20], [0  , 0  ], 'k')
ax5.text(25, -5, '$100\mu as/yr$', fontsize = 15)

plt.savefig('../plot/LD_ours.eps', dpi=100)
plt.close()