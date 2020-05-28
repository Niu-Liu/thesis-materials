# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:22:02 2015

@author: Neo
"""
#!/bin/python
#statics  file analysis

##ICRF2 defining sources data file
icrf2='../catalog/icrf2-defining.dat'
defn=[]
deff=open(icrf2)
conf=deff.readlines()
for i in range(20,len(conf)):
    dlin=str(conf[i]).split()
    defn.append(dlin[2])
deff.close()

##Some Initial Variables Defition
dat_dir='../results/'
statf='sou.sta'
#for OPA catalog
soun=[] #radio sources name
sesn=[] #length of observation sessions
obsb=[] #start epoch of observation
obss=[] #observation span
#for ICRF2 295 defining sources

##fetch data in statics file.
f_sta=open(dat_dir+statf,'r')
stat_data=f_sta.readlines()
f_sta.close()

for i in range(3,len(stat_data)):
    data=str(stat_data[i]).split()
    soun.append(data[0])
    sesn.append(int(data[1]))
    obsb.append(float(data[2]))
    obss.append(float(data[4]))
        
################################################################
# plot
################################################################        
import numpy as np
import matplotlib.pyplot as plt
#turn the data lists into arrays
sesn=np.array(sesn)
obsb=np.array(obsb)
obss=np.array(obss)

#count length of sessions and plot
#cat= ['0','20', '40', '60', '80', '100', '200', '500', '1000']
#x_pos = np.arange(0,100,10)
#y_pos = np.arange(len(cat))
#y_pos = range(len(cat)) + [len(cat)-0.5]
#y1=np.arange(0.7, len(cat), 1)
#y2=np.arange(0.3, len(cat), 1)
#
#lcou=np.zeros(len(cat),dtype=int)
#dcou=np.zeros(len(cat),dtype=int)
#
#
#for i in range(len(soun)):
##判断是否为定义源。
#    if soun[i] in defn:
#        flag=1
#    else:
#        flag=0
#        
#    num=sesn[i]
#    if num < 20:
#        lcou[0] += 1
#        dcou[0] += flag
#    elif num < 40:
#        lcou[1] += 1
#        dcou[1] += flag
#    elif num < 60:
#        lcou[2] += 1
#        dcou[2] += flag
#    elif num < 80:
#        lcou[3] += 1
#        dcou[3] += flag
#    elif num < 100:
#        lcou[4] += 1
#        dcou[4]+=flag
#    elif num < 200:
#        lcou[5]+=1
#        dcou[5]+=flag
#    elif num < 500:
#        lcou[6]+=1
#        dcou[6]+=flag
#    elif num < 1000:
#        lcou[7]+=1
#        dcou[7]+=flag   
#    else:
#        lcou[8]+=1
#        dcou[8]+=flag    
#
##what does this line mean?        
##I remember that the data file of a sources seems to be blank
##I does not know the reason.
##lcou[0] += 1
#cou1 = lcou[0]
#cou2 = lcou[1]
#lcou[0] = 200 + (lcou[0]-3000)*50.0/(3500-250)
#lcou[1] = 100 + (lcou[1]-100)/2.0

#plt.barh(y1, lcou, align='center', height=0.4, alpha=0.4, label='3826OPA')
#plt.barh(y2, dcou, align='center', height=0.4, alpha=0.4, color='y', \
#label='295ICRF')
#plt.yticks(y_pos, cat+['>1000'], fontsize=15)
#plt.xticks([0,50,100,150,200,250],('0','50','100','200','3000','3500'), fontsize=15)
#plt.ylim([0, len(cat)])
#plt.xlim([0, 250])
#plt.ylabel('Number of sources', fontsize=15)
#plt.xlabel('Number of sessions', fontsize=15)
##plt.title('Length of observation sessions')
#plt.text(lcou[0]+5, y1[0] - 0.2, str(cou1), fontsize=15)
#plt.text(dcou[0]+5, y2[0] - 0.2, str(dcou[0]), fontsize=15)
#plt.text(lcou[1]+5, y1[1] - 0.2, str(cou2), fontsize=15)
#plt.text(dcou[1]+5, y2[1] - 0.2, str(dcou[1]), fontsize=15)
#for i in range(2,len(cat)):
#    plt.text(lcou[i]+5, y1[i] - 0.2, str(lcou[i]), fontsize=15)
#    plt.text(dcou[i]+5, y2[i] - 0.2, str(dcou[i]), fontsize=15)
#plt.legend(loc='upper right', fontsize=15)
##plt.savefig('../plot/Number_of_session.ps') 
##plt.savefig('../plot/Number_of_session.eps') 
#plt.show()
##plt.close()

#Observation span 
#cat= ['0', '5', '10', '15','20','25', '30', '35']
#y_pos = range(len(cat)) + [len(cat)-0.5]
#y1=np.arange(0.7, len(cat), 1)
#y2=np.arange(0.3, len(cat), 1)
#
#lcou=np.zeros(len(cat),dtype=int)
#dcou=np.zeros(len(cat),dtype=int)
#
#for i in range(len(soun)):
##判断是否为定义源。
#    if soun[i] in defn:
#        flag=1
#    else:
#        flag=0
#        
#    spn=obss[i]
#    if spn < 5.0:
#        lcou[0]+=1
#        dcou[0]+=flag
#    elif spn < 10.0:
#        lcou[1]+=1
#        dcou[1]+=flag
#    elif spn < 15.0:
#        lcou[2]+=1
#        dcou[2]+=flag
#    elif spn < 20.0:
#        lcou[3]+=1
#        dcou[3]+=flag
#    elif spn < 25.0:
#        lcou[4]+=1
#        dcou[4]+=flag
#    elif spn < 30.0:
#        lcou[5]+=1
#        dcou[5]+=flag
#    elif spn < 35.0:
#        lcou[6]+=1
#        dcou[6]+=flag
#    else :
#        lcou[7]+=1
#        dcou[7]+=flag
#        
##lcou[0] += 1
#
#New1 = np.array(lcou)
#New2 = np.array(dcou)
#
#for i in range(len(lcou)):
#    if New1[i] > 2000:
#        New1[i] = 250 + (New1[i]-2000)*50.0/1000
#    elif New1[i] > 300:
#        New1[i] = 200 + (New1[i]-300)*50.0/1700
#    elif New1[i] > 100:
#        New1[i] = 100 + (New1[i]-100)*50.0/100
#
#for i in range(len(lcou)):
#    if New2[i] > 2000:
#        New2[i] = 250 + (New1[i]-2000)*50.0/1000
#    elif New2[i] > 300:
#        New2[i] = 200 + (New1[i]-300)*50.0/1700
#    elif New2[i] > 100:
#        New2[i] = 100 + (New1[i]-100)*50.0/100
#
#plt.barh(y1, New1, align='center', height=0.4, alpha=0.4, label='3826OPA')
#plt.barh(y2, New2, align='center', height=0.4, alpha=0.4, color='y', \
#label='295ICRF')
#plt.yticks(y_pos, cat+['>35'], fontsize=15)
#plt.xticks(np.arange(0,350,50),('0','50','100','200','300','2000','3000'), fontsize=15)
#plt.ylim([0, len(cat)])
#plt.xlim([0, 300])
#plt.xlabel('Number of sources', fontsize=15)
#plt.ylabel('Observation Span (year)', fontsize=15)
#for i in range(len(cat)):
#    plt.text(New1[i], y1[i]-0.15, str(lcou[i]), fontsize=15)
#    plt.text(New2[i], y2[i]-0.15, str(dcou[i]), fontsize=15)
#plt.legend(loc='upper right', fontsize=15)
##plt.savefig('../plot/obs-spn.png') 
#plt.savefig('../plot/Observation_span.eps')
#plt.show()
#plt.close()

#######################################################################
#######################################################################
#Cadidates Selection
#criterions: 
#1. number of onservation session is lager than or equal to 20;
#2. observation span is longer than 5 year
#
#sources should be expelled for some reasons
#* -- Radio stars
#     0236+610 0459-753 1612+339 0334+004 0323+285 1458-083 \
#* -- Known gravitational lenses
#     0218+35A 0218+35B 0218+357 1422+231 1830-21A 1830-21B 1830-211
#----39  special handling sources in IERS TN 35
#0014+813, 0106+013, 0202+149, 0208−512, 0212+735, 0235+164,
#0238−084 (NGC1052), 0316+413 (3C84), 0430+052 (3C120),
#0438−436, 0451−282, 0528+134, 0607−157, 0637−752, 0711+356,
#0738+313, 0919−260, 0923+392 (4C39.25), 0953+254 (OK290),
#1021−006, 1044+719, 1226+023 (3C273B), 1253−055 (3C279),
#1308+326, 1404+286 (OQ208), 1448+762, 1458+718 (3C309.1),
#1611+343, 1610−771, 1641+399 (3C345), 1739+522, 2121+053,
#2128−123, 2134+004, 2145+067, 2201+315, 2234+282, 2243−123,
#and 2251+158 (3C454.3).
#
#   pick up the new candidates.

exp=['0236+610', '0459-753', '1612+339', '0334+004', '0323+285', '1458-083',\
'0218+357', '1422+231', '1830-211', '0014+813', '0106+013', '0202+149',\
'0208−512', '0212+735', '0235+164', '0238−084', '0316+413', '0430+052',\
'0438−436', '0451−282', '0528+134', '0607−157', '0637−752', '0711+356',\
'0738+313', '0919−260', '0923+392', '0953+254', '1021−006', '1044+719',\
'1226+023', '1253−055', '1308+326', '1404+286', '1448+762', '1458+718',\
'1611+343', '1610−771', '1641+399', '1739+522', '2121+053', '2128−123',\
'2134+004', '2145+067', '2201+315', '2234+282', '2243−123', '2251+158']
#exp = []
cad = []

cadf = '../catalog/new_candidate1.cat'
exf =  '../catalog/ex_def1.cat'
fcad = open(cadf,'w')
fexf = open(exf, 'w')
cad = 0
for i in range(len(soun)):
#    if soun[i] in defn:
#        print >>fcad, soun[i]
#        cad += 1
#    elif sesn[i] >= 20 and obss[i] >= 10.0 and soun[i] not in exp:
#        cad += 1
#        print >>fcad, soun[i]
    if sesn[i] >= 20 and obss[i] >= 10.0 and soun[i] not in exp:
        print >>fcad, soun[i]
        cad += 1
    elif soun[i] in defn:
        print >>fexf, soun[i]

print cad     
fcad.close()
fexf.close()
