#!/usr/bin/python
## Niu, Oct 13 2016
## get the radio sources list provided on http://ivsopar.obspm.fr/radiosources/index.php

import numpy as np


s1, s2, l = np.loadtxt('index.php', dtype=str, delimiter='>', unpack=True)
#l = np.loadtxt('index.php', dtype=str, delimiter='>')
#print l[0]

fou = open('opa.list', 'w')
print>>fou, '## radio sources list provided by opa.'
for i in range(len(l)):
    print>>fou, l[i]

fou.close()

print 'Done!'
