# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 16:57:01 2016

APM / its uncertainity
let's say, Normalized Apparent proper motion


@author: neo
"""

import numpy as np
import matplotlib.pyplot as plt
from Data_Load import Apm_Load 

dat_fil = '635_sou.apm'

[Sou, V, V_err, V_RA, V_Dec, Va_E, Vd_E, RA, Dec] = Apm_Load(dat_fil)

plt.plot(V_err, V, '.')
plt.ylim([0, 200])
plt.show()