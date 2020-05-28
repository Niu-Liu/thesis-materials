# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 16:48:50 2016

plot for Deviations

@author: neo
"""

import matplotlib.pyplot as plt
import numpy as np

soun = []
Stda = []
Stdd = []
ASta = []
AStd = []
Ma   = []
Md   = [] 
WASa = []
WASd = []

N = []

in_fil = '../results/variance.dat'
fin = open(in_fil,'r')
datl = fin.readline()
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
    
    datl = fin.readline()

fin.close()    

#histogram
#y = Stdd
#spin = [0, 1, 2, 3, 4, 5, 10, 20, 50, 100]
#n, bins, pa = plt.hist(y, spin, histtype='bar', rwidth=0.8)
#plt.close()
#n = n*100.0/len(y)
#num = list(n)+ [100-sum(n)]
#x = np.arange(len(spin)) + 0.5
#xlab = [str(spin[i]) for i in range(len(spin))] + ['>'+str(spin[-1])]
#plt.bar(x, num, align='center', alpha=0.4)
#plt.xticks(np.arange(len(spin)+1), xlab)
#plt.yticks([0,10,20,30,35], ('0','10','20','30','%'))
#plt.xlabel("Weighted Standard deviation of R.A.cos(Dec)(mas)")
#plt.xlabel("Weighted Allan Standard deviation of R.A.cos(Dec)(mas)")
#plt.yticks([0,10,20,30,40,50,60,70,80], ('0','10','20','30','40','50','60','70','%'))
#plt.yticks([0,10,20,30,40,50,60], ('0','10','20','30','40','50','%'))
#plt.xlabel("Weighted Standard deviation of Declination(mas)")
#plt.yticks([0,10,20,30,40,50], ('0','10','20','30','40','%'))
#plt.xlabel("Weighted Allan Standard deviation of Declination(mas)")

#count
spin = [0, 1, 2, 5, 10, 20, 50, 100]
x = np.arange(len(spin)) + 0.5
#xlab = [str(spin[i]) for i in range(len(spin))] + ['>'+str(spin[-1])]
xtick = np.arange(len(spin))
xlab = [str(spin[i]) for i in range(len(spin))]
N = len(Stda)
#WStd for R.A.cos(Decl.)
x1 = Stda
n1, bins, pa = plt.hist(x1, spin)
n1 = n1*100.0/N
num1 = list(n1)+ [100-sum(n1)]
#WStd for Decl.
y1 = Stdd
n2, bins, pa = plt.hist(y1, spin)
n2 = n2*100.0/N
num2 = list(n2)+ [100-sum(n2)]
#WAStd for R.A.cos(Decl.)
x2 = WASa
n3, bins, pa = plt.hist(x2, spin)
n3 = n3*100.0/N
num3 = list(n3)+ [100-sum(n3)]
#WAStd for Decl.
y2 = Stdd
n4, bins, pa = plt.hist(y2, spin)
n4 = n4*100.0/N
num4 = list(n4)+ [100-sum(n4)]
plt.close()

fig, ax = plt.subplots(2, 2)
((ax1, ax2), (ax3, ax4)) = ax

ax1.bar(x, num1, align='center', alpha=0.3)
ax1.set_xticks(xtick)
ax1.set_xticklabels(xlab, fontsize = 12)
ax1.set_yticks([0,10,20,30,35])
ax1.set_yticklabels(('0','10','20','30','%'), fontsize = 15)
#ax1.set_xlabel("WDEV of R.A.cos(Dec)(mas)")
ax1.text(1, 30, 'WDEV of R.A.cos(Dec)', fontsize = 15)

ax2.bar(x, num2, align='center', alpha=0.3)
ax2.set_xticks(xtick)
ax2.set_xticklabels(xlab, fontsize = 12)
ax2.set_yticks([0,10,20,30,40,50,60])
ax2.set_yticklabels(('0','10','20','30','40','50','%'), fontsize = 15)
ax2.text(3, 50,"WDEV of Dec.",fontsize = 15)

ax3.bar(x, num3, align='center', alpha=0.3)
ax3.set_xticks(xtick)
ax3.set_xticklabels(xlab, fontsize = 12)
ax3.set_yticks([0,10,20,30,40,50,60])
ax3.set_yticklabels(('0','10','20','30','40','50','%'), fontsize = 15)
ax3.text(1, 50,"WADEV of R.A.cos(Dec)", fontsize = 15)

ax4.bar(x, num4, align='center', alpha=0.3)
ax4.set_xticks(xtick)
ax4.set_xticklabels(xlab, fontsize = 12)
ax4.set_yticks([0,10,20,30,40,50,60])
ax4.set_yticklabels(('0','10','20','30','40','50','%'), fontsize = 15)
ax4.text(3, 50,"WADEV of Dec.", fontsize = 15)

plt.savefig('../plot/DEV_plot.eps')
plt.close()