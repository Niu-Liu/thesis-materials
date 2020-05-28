# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 22:35:49 2016

@author: Neo
"""

import matplotlib.pyplot as plt
from ReadMeg import readMeg
import numpy as np
cos = np.cos
sin = np.sin

## Load proper motion data of M-K giants.
fname = '../data/merge_all_mk.fits'
#fname = '../data/merge_all_ob.fits'
tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
            plx, plx_err, l, b, g_mag = readMeg(fname, 'tycho2-id')
            
## degree -> rad.
raG  = np.deg2rad(raG)
decG = np.deg2rad(decG)
l    = np.deg2rad(l)
b    = np.deg2rad(b)
### consideration on parallax
## First eject the object with a minus plx 
indp = np.where(plx>0)[0]
## plx_err/plx ratio statistics.
rat = plx_err/plx*100     ## percentages
## set a barrier for plx_err/plx ratio
bar = 30
indr = np.where(rat<= bar)[0]
## Only the entry passing these two barriers is considered.
indplx = np.intersect1d(indp, indr)

tag    = np.take(tag,    indplx)
pmraT  = np.take(pmraT,  indplx)
pmdecT = np.take(pmdecT, indplx)
pmraG  = np.take(pmraG,  indplx)
pmdecG = np.take(pmdecG, indplx)
plx    = np.take(plx,    indplx)
raG    = np.take(raG,    indplx)
decG   = np.take(decG,   indplx)
l      = np.take(l,      indplx)
b      = np.take(b,      indplx)
## Heliocentric distance. Unit: kpc.
r = 1.0/plx
## Galactic rectangle coordinate。 Unit：kpc.
X = r*cos(l)*cos(b)
Y = r*sin(l)*cos(b)
Z = r*sin(b)

##############################################################################
## For M-K giants
## limit of r: r >= 0.3kpc
indr1 = np.where(r>=0.2)
indr2 = np.where(r<=1.0)
indr  = np.intersect1d(indr1, indr2)
#print>>flog, '## limit of r: r > 0.3kpc:\n number   %d'%len(indr[0])
# galactic height |z| <= 0.5kpc 
indz = np.where(np.abs(Z)<=0.5) 
#print>>flog, '## limit of z: |z| <= 0.5kpc:\n number   %d'%len(indz[0])
ind = np.intersect1d(indr, indz)
#print>>flog, '## Combine these two limitations:\n number   %d'%len(ind)

# extract the sub-sample.
smp = np.take(tag, ind)

rs = np.take(r, ind)
ls = np.take(l, ind)
bs = np.take(b, ind)
Xs = np.take(X, ind)
Ys = np.take(Y, ind)           
Zs = np.take(Z, ind)

indn = []

for i in range(len(tag)):
    if tag[i] not in smp:
        indn.append(i)

Xn = np.take(X, indn)
Yn = np.take(Y, indn)           
Zn = np.take(Z, indn)

# X-Y plane
plt.figure(figsize=(8,8))
A = 1.0
t = np.linspace(0, 2*np.pi, 50)
x = A*cos(t)
y = A*sin(t)
plt.plot(x,y,'g-', markersize=1)
#plt.vlines(x=0, ymin=-2.5, ymax= 2.5, linestyles='dashed')
#plt.hlines(y=0, xmin=-2.5, xmax= 2.5, linestyles='dashed')
plt.text(-0.5, -1.2, '1.0 kpc', fontsize=20)
plt.plot(Xn, Yn, 'r.', markersize=1.5)
#plt.plot(Xn, Yn, 'b.', markersize=1)
plt.plot(Xs, Ys, 'b.', markersize=1.5)
plt.xticks(fontsize=15)
plt.yticks(np.arange(-1.0, 2.0, 0.5), fontsize=15)
#A = 0.20
#x = A*cos(t)
#y = A*sin(t)
#plt.plot(x,y,'g-', markersize=1)
#A = 0.70
#x = A*cos(t)
#y = A*sin(t)
#plt.plot(x,y,'r-', markersize=1)
lim = 1.5
plt.xlim([-lim, lim])
plt.ylim([-lim, lim])
plt.xlabel('$X_G$(kpc)', fontsize=20)
plt.ylabel('$Y_G$(kpc)', fontsize=20)
plt.savefig('../plot/XY_mk.eps',dpi=500)
plt.close()

## X-Z plane
#plt.figure(figsize=(8,8))
#plt.hlines(y= 0.5, xmin=-2.5, xmax= 2.5, linestyles='dashed')
#plt.hlines(y=-0.5, xmin=-2.5, xmax= 2.5, linestyles='dashed')
#plt.text(0.9, 0.8, 'z= 0.5kpc')
#plt.text(0.9,-0.7, 'z=-0.5kpc')
#plt.plot(Xn, Zn, 'r.', markersize=1)
#plt.plot(Xs, Zs, 'b.', markersize=1)
#lim = 1.5
#plt.xlim([-lim, lim])
#plt.ylim([-lim, lim])
#plt.xlabel('X(kpc)')
#plt.ylabel('Z(kpc)')
#plt.savefig('../plot/XZ_mk.eps',dpi=500)
#plt.close()

## Y-Z plane
#plt.figure(figsize=(8,8))
#plt.hlines(y= 0.5, xmin=-2.5, xmax= 2.5, linestyles='dashed')
#plt.hlines(y=-0.5, xmin=-2.5, xmax= 2.5, linestyles='dashed')
#plt.text(0.9, 0.8, 'z= 0.5kpc')
#plt.text(0.9,-0.7, 'z=-0.5kpc')
#plt.plot(Yn, Zn, 'r.', markersize=1)
#plt.plot(Ys, Zs, 'b.', markersize=1)
#lim = 1.5
#plt.xlim([-lim, lim])
#plt.ylim([-lim, lim])
#plt.xlabel('Y(kpc)')
#plt.ylabel('Z(kpc)')
#plt.savefig('../plot/YZ_mk.eps',dpi=500)
#plt.close()

###############################################################################
## For all O-B5 data
#fname = '../data/merge_all_ob.fits'
#tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
#            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
#            plx, plx_err, l, b, g_mag = readMeg(fname, 'tycho2-id')
#            
### degree -> rad.
#raG = np.deg2rad(raG)
#decG = np.deg2rad(decG)
#l = np.deg2rad(l)
#b= np.deg2rad(b)
### Heliocentric distance. Unit: kpc.
#r = 1.0/plx
### Galactic rectangle coordinate。 Unit：kpc.
#tag1 = tag
#X = r*cos(l)*cos(b)
#Y = r*sin(l)*cos(b)
#Z = r*sin(b)
#
######################################
### For older O-B5 data
#fname = '../data/merge_ob.fits'
#tag, raT, ra_errT, decT, dec_errT, pmraT, pmra_errT, pmdecT, pmdec_errT,\
#            raG, ra_errG, decG, dec_errG, pmraG, pmra_errG, pmdecG, pmdec_errG,\
#            plx, plx_err, l, b, g_mag = readMeg(fname, 'tycho2-id')
#            
### degree -> rad.
#raG = np.deg2rad(raG)
#decG = np.deg2rad(decG)
#l = np.deg2rad(l)
#b= np.deg2rad(b)
### Heliocentric distance. Unit: kpc.
#r = 1.0/plx
### Galactic rectangle coordinate。 Unit：kpc.
#tag2 = tag
#X2 = r*cos(l)*cos(b)
#Y2 = r*sin(l)*cos(b)
#Z2 = r*sin(b)
#
#indn = []
#for i in range(len(tag1)):
#    if tag1[i] not in tag2:
#        indn.append(i)
#
#X1 = np.take(X, indn)
#Y1 = np.take(Y, indn)           
#Z1 = np.take(Z, indn)
#
################################################################################
### X-Y plane
#plt.figure(figsize=(8,8))
#m = 5
#A = 1
#t = np.linspace(0, 2*np.pi, 50)
#x = A*cos(t)
#y = A*sin(t)
#plt.plot(x,y,'g-', markersize=1)
##plt.text(1.2, 0.2, '1.0 kpc')
#plt.vlines(x=0, ymin=-m, ymax= m, linestyles='dashed')
#plt.hlines(y=0, xmin=-m, xmax= m, linestyles='dashed')
#plt.plot(X2, Y2, 'b.', markersize=1)
#plt.plot(X1, Y1, 'r.', markersize=1)
#plt.xlim([-m, m])
#plt.ylim([-m, m])
#plt.xlabel('X(kpc)')
#plt.ylabel('Y(kpc)')
#plt.savefig('../plot/XY_ob.eps',dpi=100)
#plt.close()
#
## X-Z plane
#plt.figure(figsize=(8,8))
#m = 5
#hz = 0.3
#plt.hlines(y=-hz, xmin=-m, xmax= m, linestyles='dashed')
#plt.hlines(y= hz, xmin=-m, xmax= m, linestyles='dashed')
#plt.text(-m+0.5, hz+0.2, 'z= %.1fkpc'%hz)
#plt.text(-m+0.5,-hz-0.3, 'z=-%.1fkpc'%hz)
#plt.plot(X2, Z2, 'b.', markersize=1)
#plt.plot(X1, Z1, 'r.', markersize=1)
#plt.xlim([-m, m])
#plt.ylim([-m, m])
#plt.xlabel('X(kpc)')
#plt.ylabel('Z(kpc)')
#plt.savefig('../plot/XZ_ob.eps',dpi=100)
#plt.close()
#
### Y-Z plane
#plt.figure(figsize=(8,8))
#m = 5
#hz = 0.3
#plt.hlines(y=-hz, xmin=-m, xmax= m, linestyles='dashed')
#plt.hlines(y= hz, xmin=-m, xmax= m, linestyles='dashed')
#plt.text(-m+0.5, hz+0.2, 'z= %.1fkpc'%hz)
#plt.text(-m+0.5,-hz-0.3, 'z=-%.1fkpc'%hz)
#plt.plot(Y2, Z2, 'b.', markersize=1)
#plt.plot(Y1, Z1, 'r.', markersize=1)
#plt.xlim([-m, m])
#plt.ylim([-m, m])
#plt.xlabel('Y(kpc)')
#plt.ylabel('Z(kpc)')
#plt.savefig('../plot/YZ_ob.eps',dpi=100)
#plt.close()

print 'Done!'