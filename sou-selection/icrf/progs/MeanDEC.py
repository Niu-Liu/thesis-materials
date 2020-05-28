# -*- coding: utf-8 -*-
"""
Created on Tue May 03 00:44:56 2016

@author: Neo
"""

from Data_Load import Apm_Load 
import numpy as np

#dat_fil = '394_sou.apm'
#dat_fil = '371_sou.apm'
dat_fil = '322_sou.apm'
[Sou, V, V_err, V_RA, V_Dec, EV_RA, EV_Dec, RA, Dec] = Apm_Load(dat_fil)
print np.mean(Dec)