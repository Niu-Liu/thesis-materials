# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 22:58:34 2016

@author: neo
"""

import os

#url='http://ivsopar.obspm.fr/radiosources/series/' # website of database.
url='http://ivsopar.obspm.fr/radiosources/images/'
# path to the catalog.
cat_dir='../catalog/'
#cat='sou.opa'
cat = 'OPA.cat'
dat_dir='../plot/OPA/'

f_ob = open(cat_dir+cat)
lin=f_ob.readline()

while len(lin):
    sou = lin.strip('\n')

    for k in range(len(sou)):
        if '0'<= sou[0] <='2' and sou[4]=='+':
                soun = sou[:4]+'%2B'+sou[5:]
        else:
                soun = sou
                    
    com='wget -q -N -O '+dat_dir+sou+'.png '+url+soun+'.png'
    os.system(com)
    print 'Plots for '+sou+' are Downloaded successfully.'
    
    lin=f_ob.readline()
                
f_ob.close()

print 'Download data completely.'