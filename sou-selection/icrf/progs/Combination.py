# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 21:00:13 2016

@author: Neo
"""

from Data_Load import Rot_Load 
import numpy as np
import matplotlib.pyplot as plt

#for the three special catalog
w3 = np.array([\
#icrf1    
17.81,\
#icrf2    
11.31,\
#MFV247    
10.15\
])

g3 = np.array([\
#icrf1    
0.67,\
#icrf2    
0.33,\
#MFV247    
0.69 \
])

##Here only w is needed.
[ ind, w1, wx, wy, wz ] = Rot_Load('OV.dat')
[ ind, g1, wx, wy, wz ] = Rot_Load('OVsim.dat')

[ ind, w2, wx, wy, wz ] = Rot_Load('GR.dat')
[ ind, g2, wx, wy, wz ] = Rot_Load('GRsim.dat')

w=max(max(w1),max(w2),max(w3))
g=max(max(g1),max(g2),max(g3))

t3 = w3/w*0.75 + g3/g*0.25

t1 = np.array(w1)/w*0.75 + np.array(g1)/g*0.25
t2 = np.array(w2)/w*0.75 + np.array(g2)/g*0.25

i = ind
y = np.ones(len(i))

plt.plot(i, t1, 'b')
plt.plot(i, t2, 'r')
plt.plot(i, t3[0]*y, ':', label = '212 ICRF' )
plt.plot(i, t3[1]*y, '-.', label = '295 ICRF' )
plt.plot(i, t3[2]*y, '--', label = '247 MFV' )
plt.xlim([200, max(i)])
plt.ylabel('$Q$',fontsize = 25)
plt.xlabel('No. Sources',fontsize = 15)
plt.plot(394, t2[294]-0.02, '^r', markersize=10)
plt.plot(371, t2[271]-0.02, '^r', markersize=10)
plt.plot(351, t2[251]-0.02, '^r', markersize=10)
plt.plot(322, t2[222]-0.02, '^r', markersize=10)
plt.legend()
plt.show()