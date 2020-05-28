#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 22:55:14 2016

@author: Neo

Comparison of the Galactic aberration effect between Xu, 2012 and Mignard 2016.

"""

import numpy as np
from coord_trans import Ga2Eq_rect, Eq2Ga_rect, cord2arc
import time

### The result of Xu et al 2012
### in galactic coordinates, unit: mm/s/yr
#a = np.array([ 7.47, 0.17, 3.95])
### mm/s/ur -> mas/yr
#coef = 1.0*206265e3/3.0e11
#AG = a*coef
#
### equatorial-to-Galactic matrix
#AE = Ga2Eq_rect(AG)
### accumulative effect among 13.7 year.
#g = AE*13.7
#
##print g
##g = [-0.03535305 -0.0695399  -0.01588234]

G = np.array([[ 0.00, -0.11, -0.05], \
              [-0.02, -0.08, -0.10], \
              [ 0.01, -0.08, -0.13], \
              [ 0.01, -0.07, -0.11], \
              [-0.04, -0.14, -0.07], \
              [-0.01, -0.12, -0.12]])


flog = open('../logs/GalacticAberration.log', 'a')
print>>flog, '#############################################################'
print>>flog, time.strftime('##The scripts runs at %Y-%m-%d %H:%M:%S\n',time.localtime(time.time()))

for i in range(G.size/3):
    gE = G[i]
    print>>flog, '## The Glide in equatorial coordinates (mas):\n', gE
    A = np.sqrt(np.sum(gE**2))/13.7*1000 ## uas/yr
    print>>flog, '## magnitude (uas/yr):\n', A
    gG = np.array(Eq2Ga_rect(gE)).reshape(3,)
    print>>flog, '## The Glide in galactic coordinates (mas):\n', gG
    l, b = cord2arc(gG)
    print>>flog, 'The dipole apex in galactic coordinates (l,b): (%.2f, %.2f)\n'%(np.rad2deg(l), np.rad2deg(b))