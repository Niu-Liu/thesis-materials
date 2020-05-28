'''
#!/bin/python
#Download time series data files of radio sources.
#Liu Niu, Sep 2015.
Oct 13, 2016, Niu: use numpy.loadtxt to read list file, and add the condition 
        consistent with new source names such as 'J1550-05'.
'''

import os
import numpy as np

url='http://ivsopar.obspm.fr/radiosources/series/' # website of database.

# path to the catalog.
cat_dir='../list/'
cat = 'opa.list'
dat_dir='../data/opa/'

lin = np.loadtxt(cat_dir+cat, dtype=str)

for i in range(len(lin)):
    sou = lin[i]

    if sou[0] == 'J':
        j = 5
    else:
        j = 4
        
    if sou[j]=='+':
        soun = sou[:j]+'%2B'+sou[j+1:]
    else:
        soun = sou
                    
    com='wget -q -N -O '+dat_dir+sou+'.dat '+url+soun+'.txt'
    os.system(com)
    print('Data of '+sou+' is Downloaded successfully.')
    
print('Download data completely.')