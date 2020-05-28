# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 16:48:50 2016

plot for Deviations

@author: neo

Oct 17 2016, Niu: some changes.
Oct 19 2016, Niu: change the scale of y-axis.
"""

import matplotlib.pyplot as plt
import numpy as np

in_fil = '../results/variance.dat'
soun = np.loadtxt(in_fil, dtype=str, usecols=(0,))
Stda, Asta, WASa, Stdd, AStd, WASd = \
    np.loadtxt(in_fil, usecols=list(range(3,6))+list(range(7,10)), unpack=True)   

#count
#spin = [0, 1, 2, 5, 10, 20, 50, 100]
spin = [0, 1, 2, 3, 4, 5, 10]
x = np.arange(len(spin)) + 0.5
#xlab = [str(spin[i]) for i in range(len(spin))] + ['>'+str(spin[-1])]
xtick = np.arange(len(spin)+1)
xlab = [str(spin[i]) for i in range(len(spin))]+['>10']

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
y2 = WASd
n4, bins, pa = plt.hist(y2, spin)
n4 = n4*100.0/N
num4 = list(n4)+ [100-sum(n4)]

## set y-axis scales.
yt = list(range(0, 20, 10)) + list(range(20, 50, 5)) + list(range(50, 71, 10))
ytl= ('0','10','20','','','','','','80','90','%')
## change the value of array to be consistent with the y-axis
num1[0] = num1[0]- 30.0
num2[0] = num2[0]- 30.0
num3[0] = num3[0]- 30.0
num4[0] = num4[0]- 30.0

fig, axes = plt.subplots(nrows=2, ncols=2)
ax1, ax2, ax3, ax4 = axes.flat

## now start plotting.
ax1.bar(x, num1, align='center', alpha=0.3)
ax1.set_xticks(xtick)
ax1.set_xticklabels(xlab, fontsize = 12)
ax1.set_yticks(yt)
ax1.set_yticklabels(ytl, fontsize = 12)
#ax1.set_yticks([0,10,20,30,35])
#ax1.set_yticklabels(('0','10','20','30','%'), fontsize = 15)
ax1.set_xlabel("WDEV of R.A.cos(Dec)(mas)", fontsize = 13)
ax1.set_ylabel('Percentage', fontsize = 13)
#ax1.text(1.0, 60, 'WDEV of R.A.cos(Dec)', fontsize = 15)
#ax1.axvline(x=6,hold=None,label="1",color='black',linestyle="--")
ax1.axvline(x=6,color='red',linestyle="--")

ax2.bar(x, num2, align='center', alpha=0.3)
ax2.set_xticks(xtick)
ax2.set_xticklabels(xlab, fontsize = 12)
ax2.set_yticks(yt)
ax2.set_yticklabels(ytl, fontsize = 12)
#ax2.text(1.0, 60, "WDEV of Dec.",fontsize = 15)
ax2.axvline(x=6,color='red',linestyle="--")
ax2.set_xlabel("WDEV of Dec.(mas)", fontsize = 13)
ax2.set_ylabel('Percentage', fontsize = 13)

ax3.bar(x, num3, align='center', alpha=0.3)
ax3.set_xticks(xtick)
ax3.set_xticklabels(xlab, fontsize = 12)
ax3.set_yticks(yt)
ax3.set_yticklabels(ytl, fontsize = 12)
#ax3.set_yticks([0,10,20,30,40,50,60])
#ax3.set_yticklabels(('0','10','20','30','40','50','%'), fontsize = 15)
#ax3.text(0.9, 60, "WADEV of R.A.cos(Dec)", fontsize = 15)
ax3.axvline(x=6,color='red',linestyle="--")
ax3.set_xlabel("WADEV of R.A.cos(Dec)(mas)", fontsize = 13)
ax3.set_ylabel('Percentage', fontsize = 13)

ax4.bar(x, num4, align='center', alpha=0.3)
ax4.set_xticks(xtick)
ax4.set_xticklabels(xlab, fontsize = 12)
ax4.set_yticks(yt)
ax4.set_yticklabels(ytl, fontsize = 12)
#ax4.set_yticks([0,10,20,30,40,50,60])
#ax4.set_yticklabels(('0','10','20','30','40','50','%'), fontsize = 15)
#ax4.text(1.0, 60, "WADEV of Dec.", fontsize = 15)
ax4.axvline(x=6,color='red',linestyle="--")
ax4.set_xlabel("WADEV of Dec.(mas)", fontsize = 13)
ax4.set_ylabel('Percentage', fontsize = 13)
plt.tight_layout(pad=0.5)
#plt.show()
plt.savefig('../plot/DEV_plot.eps',dpi=100)
plt.close()