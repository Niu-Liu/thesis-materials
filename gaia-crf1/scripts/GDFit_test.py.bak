# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:21:09 2016

@author: Neo

Galactic Dynamics parameter estimations.

Oct 11: using 9 parameters fitting in Subroutine ParamerFit1.
Nov  9: continue iteration until converged
"""

import numpy as np
import time
from ReadMeg import readMeg
from PMFit import PMfit, PMfit1
cos = np.cos
sin = np.sin
#import matplotlib.pyplot as plt
from coord_trans import pmT
#from OLeqnFit import ParFit

## Some constants
R0 = 8.34     ## kpc
k =4.74047

############################################################################### 
###############################################################################
def CalcRes(pml, pmb, l, b, r, x, sig):
    '''
    Calculate the residuals.
    '''
    [S1, S2, S3, A, B] = x
    [eS1, eS2, eS3, eA, eB] = sig
    
##  Calculate the residuals
    pmlc = (sin(l)*S1-cos(l)*S2)/r + A*cos(2*l)*cos(b) + B*cos(b)
    pmbc = (sin(b)*cos(l)*S1+sin(b)*sin(l)*S2-cos(b)*S3)/r -A*sin(2*l)*cos(b)*sin(b)
        
## residuals
    pmls = pml*k - pmlc
    pmbs = pmb*k - pmbc  
    varl = np.sum(pmls**2)/(len(pmls)-1)
    varb = np.sum(pmbs**2)/(len(pmbs)-1)
    wrms = np.sqrt((varl + varb)/2)
    
    return pmls, pmbs, wrms

def Elim(pmRAs,pmDEs, e_pmRA,e_pmDE, n):
    a = np.abs(pmRAs) - n*e_pmRA
    b = np.abs(pmDEs) - n*e_pmDE
    indice1 = np.where(a<0)
    indice2 = np.where(b<0)
    indice  = np.intersect1d(indice1, indice2)   

    return indice

def ParameterFit(pml, pmb, l, b, r, flog):   
## Now do the iteration until converged.    
    Num1 = 0   # number of outliers in the last iteration
    Num2 = 1   # number of outliers in this iteration
    cout = 0   # number of iterations
    n = 2.6    # n-sigma filter
    wrms = 1   # mean uncertainty
    
    npml, npmb, nl, nb, nr = pml, pmb, l, b, r
        
    while Num1 != Num2:
        Num1 = Num2
        x, sig, corrcoef = \
            PMfit(npml,npmb,wrms*np.ones(len(npml)),wrms*np.ones(len(npmb)),nl,nb,nr)
        pmls, pmbs, wrms  = CalcRes(pml, pmb, l, b, r, x, sig) 
## generate the new array after eliminating outliers 
        indice = Elim(pmls,pmbs, wrms*np.ones(len(pmls)),wrms*np.ones(len(pmbs)), n)
        npml = np.take(pml, indice)
        npmb = np.take(pmb, indice)
        nl   = np.take(l,   indice)
        nb   = np.take(b,   indice) 
        nr   = np.take(r,   indice)
        Num2 = len(pml)-len(npml)
        cout += 1
        print>>flog, 'Iteration %d: %d'%(cout, Num2)
        
    print>>flog, 'Eliminating(%2.1f-sigma): %d/%d (objects/Outliers).'%(n,len(pml), Num2)
## final result.
    '''
    (S1, S2, S3) -> peculiar motion of Sun relative the LSR
    A,B -> Oort constants
    '''    
    [S1, S2, S3, A, B] = x
    [eS1, eS2, eS3, eA, eB] = sig
        
    print>>flog, '''
    After elimination: '''
    print>>flog, 'The parameters S1, S2, S3(km/s):\n', \
        '    %.2f \pm %.2f'*3%(S1, eS1, S2, eS2, S3, eS3)
    print>>flog, ' Oort constants A, B(km/s/kpc) :\n ',  \
        '    %.2f \pm %.2f'*2%(A, eA, B, eB)
    print>>flog, ' Oort constants A, B(mas/yr) :\n ',  \
        '    %.3f \pm %.3f'*2%(A/k, eA/k, B/k, eB/k)
    print>>flog, '   correlation coefficients are:\n', corrcoef
    
    V0 = (A-B)*R0
    eV0 = R0*np.sqrt(eA**2 + eB**2)
    print>>flog, 'rotational velocity: %.2f \pm %.2f km/s.\n'%(V0, eV0)
############################################################################### 
###############################################################################  
def CalcRes1(pml, pmb, l, b, r, x, sig):
    '''
    Calculate the residuals.
    '''
    [S1, S2, S3, D32M, D13M, D21M, D12P, D13P, D32P] = x
    [eS1, eS2, eS3, eD32M, eD13M, eD21M, eD12P, eD13P, eD32P] = sig
    
##  Calculate the residuals        
    pmlc = (sin(l)*S1 - cos(l)*S2 + 0)/r + \
            D32M*(-cos(l)*sin(b)) + D13M*(-sin(l)*sin(b)) + D21M*cos(b) + \
            D12P*(cos(2*l)*cos(b)) + D13P*(-sin(l)*sin(b)) + D32P*(cos(l)*sin(b))
    pmbc = (sin(b)*cos(l)*S1+sin(b)*sin(l)*S2-cos(b)*S3)/r + \
            D32M*sin(l) + D13M*(-cos(l)) + 0.0 + \
            D12P*(-sin(2*l)*sin(b)*cos(b)) + D13P*(cos(l)*cos(2*b)) + D32P*(sin(l)*cos(2*b))
    
## residuals
    pmls = pml*k - pmlc
    pmbs = pmb*k - pmbc  
    varl = np.sum(pmls**2)/(len(pmls)-1)
    varb = np.sum(pmbs**2)/(len(pmbs)-1)
    wrms = np.sqrt((varl + varb)/2)
    
    return pmls,pmbs, wrms
    
def ParameterFit1(pml, pmb, l, b, r, flog):
## Now do the iteration until converged. 
    Num1 = 0   # number of outliers in the last iteration
    Num2 = 1   # number of outliers in this iteration
    cout = 0   # number of iterations
    n = 2.6    # n-sigma filter
    wrms = 1   # mean uncertainty
    
    npml, npmb, nl, nb, nr = pml, pmb, l, b, r
        
    while Num1 != Num2:
        Num1 = Num2
        x, sig, corrcoef = \
            PMfit1(npml, npmb, wrms*np.ones(len(npml)), wrms*np.ones(len(npmb)),nl,nb,nr)
        pmls,pmbs, wrms  = CalcRes1(pml, pmb, l, b, r, x, sig)    
## generate the new array after eliminating outliers 
        indice = Elim(pmls, pmbs, wrms*np.ones(len(pmls)), wrms*np.ones(len(pmbs)), n)
        npml = np.take(pml, indice)
        npmb = np.take(pmb, indice)
        nl = np.take(l, indice)
        nb = np.take(b, indice) 
        nr = np.take(r, indice)
        Num2 = len(pml)-len(npml)
        cout += 1
        print>>flog, 'Iteration %d: %d'%(cout, Num2)
        
    print>>flog, 'Eliminating(%2.1f-sigma): %d/%d (objects/Outliers).'%(n, len(pml), Num2)
## final result. 
    '''
    (S1, S2, S3) -> peculiar motion of Sun relative the LSR
    (D32M, D13M, D21M) -> vortivity
    (D12P, D13P, D32P) -> shear
    '''    
    [S1, S2, S3, D32M, D13M, D21M, D12P, D13P, D32P] = x
    [eS1, eS2, eS3, eD32M, eD13M, eD21M, eD12P, eD13P, eD32P] = sig
    print>>flog, 'The parameters S1, S2, S3(km/s):\n', \
        '    %.2f \pm %.2f'*3%(S1, eS1, S2, eS2, S3, eS3)
    print>>flog, ' Rotation D32M, D13M, D21M (km/s/kpc) :\n ',  \
        '    %.2f \pm %.2f'*3%(D32M, eD32M, D13M, eD13M, D21M, eD21M)
    print>>flog, ' Shear D12P, D13P, D32P (km/s/kpc) :\n ',  \
        '    %.2f \pm %.2f'*3%(D12P, eD12P, D13P, eD13P, D32P, eD32P)   
    print>>flog, ' Rotation D32M, D13M, D21M (mas/yr) :\n ',  \
        '    %.3f \pm %.3f'*3%(D32M/k, eD32M/k, D13M/k, eD13M/k, D21M/k, eD21M/k)
    print>>flog, ' Shear D12P, D13P, D32P (mas/yr) :\n ',  \
        '    %.3f \pm %.3f'*3%(D12P/k, eD12P/k, D13P/k, eD13P/k, D32P/k, eD32P/k)    
    print>>flog, '   correlation coefficients are:\n', corrcoef
            
    V0 = (D12P - D21M)*R0
    eV0 = R0*np.sqrt(eD12P**2 + eD21M**2)
    
    print>>flog, 'rotational velocity: %.2f \pm %.2f km/s.\n'%(V0, eV0)
    
############################################################################### 
###  Main Body
###############################################################################
## Log file.
flog = open('../logs/GDfit.log', 'a')
print>>flog, time.strftime('\nThe scripts runs at %Y-%m-%d %H:%M:%S',time.localtime(time.time()))

## Load proper motion data of M-K giants.
fname = '../data/merge_all_mk.fits'
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
bar = 100
indb = np.where(rat<= bar)[0]
#print indr
## Only the entry passing these two barriers is considered.
indplx = np.intersect1d(indp, indb)

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

## limit of r: r >= 0.3kpc
indr1 = np.where(r>=0.3)
indr2 = np.where(r<=1.0)
indr  = np.intersect1d(indr1, indr2)
print>>flog, '## limit of r: 1.0 > r > 0.3kpc:\n number   %d'%len(indr)
## galactic height |z| <= 0.5kpc 
indz = np.where(np.abs(Z)<=0.5) 
print>>flog, '## limit of z: |z| <= 0.5kpc:\n number   %d'%len(indz[0])
ind = np.intersect1d(indr, indz[0])
print>>flog, '## Combine these two limitations:\n number   %d'%len(ind)

## extract the sub-sample.
tags    = np.take(tag,    ind)
pmraTs  = np.take(pmraT,  ind)
pmdecTs = np.take(pmdecT, ind)
pmraGs  = np.take(pmraG,  ind)
pmdecGs = np.take(pmdecG, ind)
raGs    = np.take(raG,    ind)
decGs   = np.take(decG,   ind)
plxs    = np.take(plx,    ind)
rs = np.take(r, ind)
ls = np.take(l, ind)
bs = np.take(b, ind)
Xs = np.take(X, ind)
Ys = np.take(Y, ind)           
Zs = np.take(Z, ind)
#rMs = 1.0/(plxs-0.3)
#rPs = 1.0/(plxs+0.3)

## proper motion from 
pmlTs, pmbTs, e_pmlTs, e_pmbTs = pmT(pmraTs, pmdecTs, np.ones(len(pmraTs)), np.ones(len(pmraTs)),\
raGs, decGs, ls)
pmlGs, pmbGs, e_pmlGs, e_pmbGs = pmT(pmraGs, pmdecGs, np.ones(len(pmraGs)), np.ones(len(pmraGs)),\
raGs, decGs, ls)

###Parameter Fitting.
#print>>flog, '''\n#############################################################
### Using Tycho-2 catalog:'''
#print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
#ParameterFit(pmlTs, pmbTs, ls, bs, rs, flog)
#ParameterFit1(pmlTs, pmbTs, ls, bs, rs, flog)
#ParameterFit2(pmlTs, ls, bs, rs, flog)

print>>flog, '''\n#############################################################
## Using TGAS catalog:'''
print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
ParameterFit( pmlGs, pmbGs, ls, bs, rs, flog)
ParameterFit1(pmlGs, pmbGs, ls, bs, rs, flog)
#ParameterFit2(pmlGs, ls, bs, rs, flog)

### the systematic error is not considered temporarily, advised by Prof. Zhu. 
#print>>flog, '''\n################## consider the influence of systematic error 0.3mas
#'''
### plx = plx + 0.3
#print>>flog, '1. plx = plx + 0.3:'
#ParameterFit(pmlGs, pmbGs, ls, bs, rPs, flog)
#ParameterFit1(pmlGs, pmbGs, ls, bs, rPs, flog)
#ParameterFit2(pmlGs, ls, bs, rPs, flog)
#print>>flog, '2. plx = plx - 0.3:'
### plx = plx - 0.3
#ParameterFit(pmlGs, pmbGs, ls, bs, rMs, flog)
#ParameterFit1(pmlGs, pmbGs, ls, bs, rMs, flog)
#ParameterFit2(pmlGs, ls, bs, rMs, flog)
#
#### difference of proper motion
#dpml = pmlGs - pmlTs
#dpmb = pmbGs - pmbTs
#print>>flog, '''\n#############################################################
### Proper motion differences bewteen two catalogs:'''
#print>>flog, time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
#ParameterFit(dpml, dpmb, ls, bs, rs, flog)
#ParameterFit1(dpml, dpmb, ls, bs, rs, flog)
#ParameterFit2(dpml, ls, bs, rs, flog)

###############################################################################
### Plot for pml
#pi = np.pi
#plt.figure(figsize=(20,10))
#plt.plot(ls/pi, dpml/cos(bs), 'c.', markersize=1.5)
#plt.ylim([-20, 20])
#plt.yticks(fontsize=20)
#plt.ylabel('$\mu _l(mas/yr)$', fontsize=25, rotation=90)
#plt.xticks([0, 0.5, 1, 1.5, 2], ['0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'], \
#    fontsize=22, rotation=10)
#plt.xlabel('$Galactic\, Longitude$', fontsize=25)
#plt.hlines(y=0, xmin=0, xmax=2, linestyles='dashed', colors='r')
#A = -0.452/k
#B = -0.314/k
#x = np.linspace(0, 2, 100)
#y = A*cos(2*pi*x) + B
#plt.plot(x, y, 'r-')
#plt.hlines(y=B, xmin=0, xmax=2, linestyles='dashed', colors='k')
#plt.savefig('../plot/pml.eps', dpi=100)

print'Done!'