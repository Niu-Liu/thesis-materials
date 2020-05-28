# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 21:56:00 2016

@author: neo

Linear drift plot. 

"""

import matplotlib.pyplot as plt
from Data_Load import Apm_Load 

#dat_fil = 'icrf1_sou.apm'
#dat_fil = 'icrf2_sou.apm'
#dat_fil = 'MFV247_sou.apm'
dat_fil = 'AMS260_sou.apm'

[Sou, V, V_err, V_RA, V_Dec, Va_E, Vd_E, RA, Dec] = Apm_Load(dat_fil)
#
#dat_fil = 'Sou307_R1.dat'
#dat_fil = 'Sou366_R1.dat'
#dat_fil = 'Sou304_R2.dat'
#dat_fil = 'Sou324_R2.dat'

icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

#dat_dir =  '../results/'
#fdat = open(dat_dir + dat_fil, 'r')
#data = fdat.readlines()
#
#Sou = []    
#V = []
#V_RA  = []
#V_Dec  = []
#    
#V_err = []
#Va_E  = []
#Vd_E  = []
#    
#RA = []
#Dec = []
#
#for i in range(len(data)):
#    datl = data[i].strip('\n').split()
#    Sou.append(datl[0])
#    V.append(float(datl[1]))
#    V_RA.append(float(datl[3]))
#    V_Dec.append(float(datl[4]))
#
#    V_err.append(float(datl[2]))
#    Va_E.append(float(datl[5]))
#    Vd_E.append(float(datl[6]))
#    
#    RA.append(float(datl[7]))
#    Dec.append(float(datl[8]))
#    
print sum(Dec)/len(Dec)
        
for i in range(len(Sou)):
    if V[i] >= 500:
        plt.arrow(RA[i], Dec[i], \
        V_RA[i]/10, V_Dec[i]/10, width = 0.05, color = 'r')
    else:
        plt.arrow(RA[i], Dec[i], \
        V_RA[i]/10, V_Dec[i]/10, width = 0.05, color = 'b')
    if Sou[i] in defn:
        plt.plot(RA[i], Dec[i], 'k.')
    else:
        plt.plot(RA[i], Dec[i], 'r.')

plt.xlim(0, 360)
plt.ylim(-90, 90)
plt.xticks([0, 45, 90, 135, 180, 225, 270, 315, 360],\
['0', '', '90', '', '180', '', '270', '', '360'])
plt.yticks([-90, -60, -30, 0, 30, 60, 90])

plt.xlabel('Right Ascension($\degree$)', fontsize = 20)
plt.ylabel('  Declination  ($\degree$)', fontsize = 20)
plt.grid(True)

## Proper motion scale
plt.plot([10, 10], [-80.5, -79.5], 'k')
plt.plot([20, 20], [-80.5, -79.5], 'k')
plt.plot([10, 20], [-80  , -80  ], 'k')
plt.text(8, -78, '$100\mu as/yr$', fontsize = 15)

plt.legend()
plt.show()        
