# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 21:12:01 2016

@author: neo

Generate the catalogs of defining
"""

from Data_Load import Rank_Load

#for rank1
#list1, number = 307
#num = 307
#res = 'Sou307_R1.dat'
num = 366
res = 'Sou366_R1.dat'

#for rank2
#num = 304
#res = 'Sou304_R2.dat'
#num = 324
#res = 'Sou324_R2.dat'

dat_fil = 'NC3.rank'
res_dir = '../results/'
fres = open(res_dir + res, 'w')

[Sou, V, V_err, V_RA,V_Dec, Va_E, Vd_E, RA, Dec] = Rank_Load(dat_fil)

#for i in range(num):
#    print >> fres, Sou[i] + "  %10.4f    %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f" \
#    %(V[i],V_err[i], V_RA[i], V_Dec[i], Va_E[i], Vd_E[i], RA[i], Dec[i])
l = []
z = 75
nodes = [-90, -30, 0, 30, 90]
while num != len(l):
    l = []
    for j in range(len(nodes)-1):
        a = []
        for i in range(len(V)):
            if nodes[j] < Dec[i] < nodes[j+1]:
                a.append(Sou[i])
                
        a = a[:z]
        l += a
    print len(l)
        
    z += 1
        
for i in range(len(V)):
    if Sou[i] in l:    
        print >> fres, Sou[i] +\
"  %10.4f    %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f" \
    %(V[i],V_err[i], V_RA[i], V_Dec[i], Va_E[i], Vd_E[i], RA[i], Dec[i])
    
fres.close()
print 'Done!'
