#!/bin/python
#Download time series data files of radio sources.
#Liu Niu, Sep 2015.

import os

url='http://ivsopar.obspm.fr/radiosources/series/' # website of database.

# path to the catalog.
cat_dir='../catalog/'
#cat='sou.opa'
cat = 'OPA.cat'
dat_dir='../data/opa/'

f_ob = open(cat_dir+cat)
lin=f_ob.readline()

while len(lin):
    sou = lin.strip('\n')

    for k in range(len(sou)):
        if '0'<= sou[0] <='2' and sou[4]=='+':
                soun = sou[:4]+'%2B'+sou[5:]
        else:
                soun = sou
                    
    com='wget -q -N -O '+dat_dir+sou+'.dat '+url+soun+'.txt'
    os.system(com)
    print 'Data of '+sou+' is Downloaded successfully.'
    
    lin=f_ob.readline()
                
f_ob.close()

print 'Download data completely.'
