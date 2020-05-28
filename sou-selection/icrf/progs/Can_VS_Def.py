# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:09:17 2016

@author: Neo
"""

import fun
from Data_Load import Apm_Load 

dat_fil = '718_sou.apm'
[Sou, V, V_RA, V_De, V_err, V_eRA, V_eDe, RA, De] = Apm_Load(dat_fil)

##ICRF2 defining sources data file
icrf2 = '../catalog/icrf2-defining.dat'
defn  = []
deff  = open(icrf2)
conf  = deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

cat_dir = '../catalog/'
cand    = '300can.list'
cadn = []
cadf = open(cat_dir +cand, 'r')
cadl = cadf.readlines()
cadf.close()
for i in range(len(cadl)):
    cadn.append(cadl[i].strip('\n'))
    
V_RAd = []
V_Ded = []
V_eRAd = []
V_eDed = [] 
RAd  = []
Ded = []

V_RAc = []
V_Dec = []
V_eRAc = []
V_eDec = [] 
RAc  = []
Dec = []

for i in range(len(Sou)):
    if Sou[i] in defn:
        V_RAd.append(V_RA[i])
        V_Ded.append(V_De[i])
        V_eRAd.append(V_eRA[i])
        V_eDed.append(V_eDe[i])
        RAd.append(RA[i])
        Ded.append(De[i])
        
    if Sou[i] in cadn:
        V_RAc.append(V_RA[i])
        V_Dec.append(V_De[i])
        V_eRAc.append(V_eRA[i])
        V_eDec.append(V_eDe[i])
        RAc.append(RA[i])
        Dec.append(De[i])
        
out_dir = '../results/'
out_fil = 'Axis_Rotation.txt'
fout = open(out_dir + out_fil, 'w')
print >>fout, '#w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz'
[w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi] \
= fun.axis_rot(V_RAd, V_Ded, V_eRAd, V_eDed, RAd, Ded)
print >>fout, \
"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f" \
    %(w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz)

del w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz

[w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi] \
= fun.axis_rot(V_RAc, V_Dec, V_eRAc, V_eDec, RAc, Dec)
print >>fout, \
"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f" \
    %(w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz)
    
del w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz
#print fun.mean_dec(defn, Ded)
#print fun.mean_dec(cadn, Dec)

# Let take glide into consideration
print >>fout, '##w,sigw,wx,wy,wz,sigwx,sigwy,sigwz '
print >>fout, '##g,sigg,gx,gy,gz,siggx,siggy,siggz'
[w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,g,sigg,gx,gy,gz,siggx,siggy,siggz,chi] \
= fun.Rotation_Glide(V_RAd, V_Ded, V_eRAd, V_eDed, RAd, Ded)
print >>fout, \
"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f" \
%(w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)
print >>fout, \
"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f" \
%(g,sigg,gx,gy,gz,siggx,siggy,siggz)

del w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,g,sigg,gx,gy,gz,siggx,siggy,siggz,chi

[w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,g,sigg,gx,gy,gz,siggx,siggy,siggz,chi] \
= fun.Rotation_Glide(V_RAc, V_Dec, V_eRAc, V_eDec, RAc, Dec)
print >>fout, \
"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f" \
%(w,sigw,wx,wy,wz,sigwx,sigwy,sigwz)
print >>fout, \
"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f" \
%(g,sigg,gx,gy,gz,siggx,siggy,siggz)

del w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,g,sigg,gx,gy,gz,siggx,siggy,siggz,chi

fout.close()
