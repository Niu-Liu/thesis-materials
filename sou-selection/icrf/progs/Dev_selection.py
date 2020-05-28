# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:22:54 2016

pre-selection(criterion 3)

@author: neo
"""
##ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

soun = []
Stda = []
Stdd = []
ASta = []
AStd = []
Ma   = []
Md   = [] 
WASa = []
WASd = []

N    = []
gap  = []

in_fil = '../results/variance.dat'
fin = open(in_fil,'r')
datl = fin.readline()
out_fil = '../catalog/new_candidate2.cat'
fou = open(out_fil,'w')
while len(datl):
    dat = datl.strip('\n').split()
    soun.append(dat[0])
    Ma.append(float(dat[3]))
    Stda.append(float(dat[4]))
    ASta.append(float(dat[5]))
    WASa.append(float(dat[6]))
    Md.append(float(dat[7]))
    Stdd.append(float(dat[8]))
    AStd.append(float(dat[9]))
    WASd.append(float(dat[10]))
    
    N.append(int(dat[1]))
    gap.append(int(dat[11]))
    
    datl = fin.readline()

ex = []
for i in range(len(soun)):
#    if gap[i] <= 5 or N[i] >= 12: 
    if max(Stda[i],Stdd[i],WASa[i],WASd[i]) < 100:
        print >> fou, soun[i]
    else:
        ex.append(soun[i])
        
fou.close()
    
#ex = [soun[i] for i in range(len(soun)) if gap[i] > 5 ]
#print len(ex)
defx = [ex[i] for i in range(len(ex)) if ex[i] in defn]
ex_fil = '../catalog/ex_def2.cat'
fex = open(ex_fil, 'w')
for i in range(len(defx)):
    print>>fex, defx[i]
fex.close()