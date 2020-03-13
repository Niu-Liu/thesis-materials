#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:01:56 2016

@author: Neo

different samples of M-K giants of heliocentric distance r has different lower limits.
"""

import numpy as np
import time
from ReadMeg import readMeg
from PMFit import PMfit, PMfit1, PMfit2
cos = np.cos
sin = np.sin
from coord_trans import pmT
#from OLeqnFit import ParFit

## Some constants
R0 = 8.34     ## kpc
k  = 4.74047

############################################################################### 
###############################################################################
def CalcRes(pml, pmb, l, b, r, x, sig):
    '''
    Calculate the residuals.
    '''
    [ S1,  S2,  S3,  A,  B] = x
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

def ParameterFit(pml, pmb, l, b, r, flog, flog1):   
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
## used in latex.    
    print>>flog1, '&%d/%d'%(len(pml), Num2) + \
        '\t&$%.2f \pm %.2f$'*3%(S1, eS1, S2, eS2, S3, eS3) + \
        '\t&$%.2f \pm %.2f$'*2%(A, eA, B, eB) + \
        '\t&$%.2f \pm %.2f$'%(V0, eV0)
    
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
    
#    wrmsl = np.sqrt(varl)
#    wrmsb = np.sqrt(varb)
#    
#    return pmls,pmbs, wrmsl, wrmsb
    
def ParameterFit1(pml, pmb, l, b, r, flog, flog1):
## Now do the iteration until converged. 
    Num1 = 0   # number of outliers in the last iteration
    Num2 = 1   # number of outliers in this iteration
    cout = 0   # number of iterations
    n    = 2.6 # n-sigma filter
#    wrmsl= 1.0 # mean uncertainty
#    wrmsb= 1.0
    wrms = 1.0
    
    npml, npmb, nl, nb, nr = pml, pmb, l, b, r
        
    while Num1 != Num2:
        Num1 = Num2
        x, sig, corrcoef = \
            PMfit1(npml, npmb, wrms*np.ones(len(npml)), wrms*np.ones(len(npmb)),nl,nb,nr)
        pmls,pmbs, wrms = CalcRes1(pml, pmb, l, b, r, x, sig)
#        x, sig, corrcoef = \
#            PMfit1(npml, npmb, wrmsl*np.ones(len(npml)), wrmsb*np.ones(len(npmb)),nl,nb,nr)
#        pmls,pmbs, wrmsl, wrmsb = CalcRes1(pml, pmb, l, b, r, x, sig)    
## generate the new array after eliminating outliers 
#        indice = Elim(pmls, pmbs, wrmsl*np.ones(len(pmls)), wrmsb*np.ones(len(pmbs)), n)
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
## used in latex.    
#    print>>flog1, '&%d/%d'%(len(pml), Num2) + \
#        '\t&$%.2f \pm %.2f$'*3%(S1, eS1, S2, eS2, S3, eS3) + \
#        '\t&$%.2f \pm %.2f$'*3%(D32M, eD32M, D13M, eD13M, D21M, eD21M) + \
#        '\t&$%.2f \pm %.2f$'*3%(D12P, eD12P, D13P, eD13P, D32P, eD32P) + \
#        '\t&$%.2f \pm %.2f$'%(V0, eV0) 
     
    print>>flog1, '&%d/%d'%(len(pml), Num2) + \
        '##val'+'  \t&$%.2f $'*10%(S1, S2, S3, D32M, D13M, D21M, D12P, D13P, D32P, V0) 
    print>>flog1, \
        '##err'+'  \t&$\pm %.2f$'*10%(eS1, eS2, eS3, eD32M, eD13M, eD21M,\
                                      eD12P, eD13P, eD32P, eV0) 
                
############################################################################### 
###############################################################################  
def CalcRes2(pml, pmb, l, b, r, x, sig):
    '''
    Calculate the residuals.
    '''
    [S1, S2, S3, D13M, D21M, D12P, D13P] = x
    [eS1, eS2, eS3, eD13M, eD21M, eD12P, eD13P] = sig
    
##  Calculate the residuals        
    pmlc = (sin(l)*S1 - cos(l)*S2 + 0)/r + \
             + D13M*(-sin(l)*sin(b)) + D21M*cos(b) + \
            D12P*(cos(2*l)*cos(b)) + D13P*(-sin(l)*sin(b))
    pmbc = (sin(b)*cos(l)*S1+sin(b)*sin(l)*S2-cos(b)*S3)/r + \
             + D13M*(-cos(l)) + 0.0 + \
            D12P*(-sin(2*l)*sin(b)*cos(b)) + D13P*(cos(l)*cos(2*b))
    
## residuals
    pmls = pml*k - pmlc
    pmbs = pmb*k - pmbc  
    varl = np.sum(pmls**2)/(len(pmls)-1)
    varb = np.sum(pmbs**2)/(len(pmbs)-1)
    
    wrms = np.sqrt((varl + varb)/2)
    return pmls,pmbs, wrms
    
#    wrmsl = np.sqrt(varl)
#    wrmsb = np.sqrt(varb)
#    
#    return pmls,pmbs, wrmsl, wrmsb
    
def ParameterFit2(pml, pmb, l, b, r, flog, flog1):
## Now do the iteration until converged. 
    Num1 = 0   # number of outliers in the last iteration
    Num2 = 1   # number of outliers in this iteration
    cout = 0   # number of iterations
    n    = 2.6 # n-sigma filter
#    wrmsl= 1.0 # mean uncertainty
#    wrmsb= 1.0
    wrms = 1.0
    
    npml, npmb, nl, nb, nr = pml, pmb, l, b, r
        
    while Num1 != Num2:
        Num1 = Num2
        x, sig, corrcoef = \
            PMfit2(npml, npmb, wrms*np.ones(len(npml)), wrms*np.ones(len(npmb)),nl,nb,nr)
        pmls,pmbs, wrms = CalcRes2(pml, pmb, l, b, r, x, sig)
#        x, sig, corrcoef = \
#            PMfit1(npml, npmb, wrmsl*np.ones(len(npml)), wrmsb*np.ones(len(npmb)),nl,nb,nr)
#        pmls,pmbs, wrmsl, wrmsb = CalcRes1(pml, pmb, l, b, r, x, sig)    
## generate the new array after eliminating outliers 
#        indice = Elim(pmls, pmbs, wrmsl*np.ones(len(pmls)), wrmsb*np.ones(len(pmbs)), n)
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
    (0,    D13M, D21M) -> vortivity
    (D12P, D13P,    0) -> shear
    '''    
    [S1,   S2,  S3,  D13M,  D21M,  D12P,  D13P] = x
    [eS1, eS2, eS3, eD13M, eD21M, eD12P, eD13P] = sig
    print>>flog, 'The parameters S1, S2, S3(km/s):\n', \
        '    %.2f \pm %.2f'*3%(S1, eS1, S2, eS2, S3, eS3)
    print>>flog, ' Rotation D32M, D13M, D21M (km/s/kpc) :\n ',  \
        '    %.2f \pm %.2f'*3%(0, 0, D13M, eD13M, D21M, eD21M)
    print>>flog, ' Shear D12P, D13P, D32P (km/s/kpc) :\n ',  \
        '    %.2f \pm %.2f'*3%(D12P, eD12P, D13P, eD13P, 0, 0)   
    print>>flog, ' Rotation D32M, D13M, D21M (mas/yr) :\n ',  \
        '    %.2f \pm %.2f'*3%(0/k, 0/k, D13M/k, eD13M/k, D21M/k, eD21M/k)
    print>>flog, ' Shear D12P, D13P, D32P (mas/yr) :\n ',  \
        '    %.2f \pm %.2f'*3%(D12P/k, eD12P/k, D13P/k, eD13P/k, 0/k, 0/k)    
    print>>flog, '   correlation coefficients are:\n', corrcoef
            
    V0 = (D12P - D21M)*R0
    eV0 = R0*np.sqrt(eD12P**2 + eD21M**2)
    
    print>>flog, 'rotational velocity: %.2f \pm %.2f km/s.\n'%(V0, eV0)
## used in latex.    
#    print>>flog1, '&%d/%d'%(len(pml), Num2) + \
#        '\t&$%.2f \pm %.2f$'*3%(S1, eS1, S2, eS2, S3, eS3) + \
#        '\t&$%.2f \pm %.2f$'*2%(D13M, eD13M, D21M, eD21M) + \
#        '\t&$%.2f \pm %.2f$'*2%(D12P, eD12P, D13P, eD13P) + \
#        '\t&$%.2f \pm %.2f$'%(V0, eV0)
    print>>flog1, '&%d/%d'%(len(pml), Num2) + \
        '##val'+'  \t&$%.2f $'*8%(S1, S2, S3, D13M, D21M, D12P, D13P, V0) 
    print>>flog1, \
        '##err'+'  \t&$\pm %.2f$'*8%(eS1, eS2, eS3, eD13M, eD21M,\
                                      eD12P, eD13P, eV0) 
          
################################################################################ 
################################################################################    
#def ParameterFit2(pml, l, b, r, flog):
####  First Loop, calculate the mean error without weights. 
#    print>>flog, '''###########################################################
#    ### Using only proper motion of longtitude components: '''
#   
#    x, sig, corrcoef = ParFit(pml,np.ones(len(pml)),l,b,r)
#    '''
#    (S1, S2, S3) -> peculiar motion of Sun relative the LSR
#    A,B -> Oort constants
#    '''    
#    [S1, S2, A, B] = x
#    [eS1, eS2, eA, eB] = sig
###  Calculate the residuals    
#    pmlc = (sin(l)*S1-cos(l)*S2)/r + A*cos(2*l)*cos(b) + B*cos(b)
#        
### residuals
#    pmls = pml*k - pmlc 
#    wrmsl = np.sqrt(np.sum(pmls**2)/(len(pmls)-1))
#    wrms = wrmsl
#    
### generate the new array after eliminating outliers 
#    n = 2.6
#    indice = Elim1(pmls, wrms*np.ones(len(pmls)), n)
#    npml = np.take(pml, indice)
#    nl = np.take(l, indice)
#    nb = np.take(b, indice) 
#    nr = np.take(r, indice)
#    
#    x, sig, corrcoef = ParFit(npml,np.ones(len(npml[0])),nl,nb,nr)
###  Calculate the residuals
#    [S1, S2, A, B] = x
#    [eS1, eS2, eA, eB] = sig
#    pmlc = (sin(l)*S1-cos(l)*S2)/r + A*cos(2*l)*cos(b) + B*cos(b)
#       
### residuals
#    pmls = pml*k - pmlc
#    wrmsl = np.sqrt(np.sum(pmls**2)/(len(pmls)-1))
#    wrms = wrmsl
#    merr = wrms 
#        
#    ##  Next, apply the weigths of 1/M_err**2
#    print>>flog, 'M_err = %f mas*yr-1'%merr
#    x, sig, corrcoef = ParFit(pml,merr*np.ones(len(pml)),l,b,r)
#    [S1, S2, A, B] = x
#    [eS1, eS2, eA, eB] = sig
###  Calculate the residuals
#    pmlc = (sin(l)*S1-cos(l)*S2)/r + A*cos(2*l)*cos(b) + B*cos(b)
#    
### residuals
#    pmls = pml*k - pmlc
#    
### generate the new array after eliminating outliers 
#    n = 2.6
#    indice = Elim1(pmls, merr*np.ones(len(pmls)), n)
#    npml = np.take(pml, indice)
#    nl = np.take(l, indice)
#    nb = np.take(b, indice) 
#    nr = np.take(r, indice)
#        
#    print>>flog, 'Eliminating(%2.1f-sigma): %d/%d (objects/Outliers).'%(n,len(pml), len(pml)-len(npml[0]))
#        
#    x, sig, corrcoef = ParFit(npml,merr*np.ones(len(npml[0])),nl,nb,nr)
#    [S1, S2, A, B] = x
#    [eS1, eS2, eA, eB] = sig
#    print>>flog, '''
#    After elimination: '''
#    print>>flog, 'The parameters S1, S2, S3(km/s):\n', \
#        '    %.2f \pm %.2f'*2%(S1, eS1, S2, eS2)
#    print>>flog, ' Oort constants A, B(km/s/kpc) :\n ',  \
#        '    %.2f \pm %.2f'*2%(A, eA, B, eB)  
#    print>>flog, ' Oort constants A, B(mas/yr) :\n ',  \
#        '    %.3f \pm %.3f'*2%(A/k, eA/k, B/k, eB/k) 
#    print>>flog, '   correlation coefficients are:\n', corrcoef
#    
#    V0 = (A-B)*R0
#    eV0 = R0*np.sqrt(eA**2 + eB**2)
#    print>>flog, 'rotational velocity: %.2f \pm %.2f km/s.\n'%(V0, eV0)
 
## Log file.
flog = open('../logs/DiffSamp.log', 'a')
print>>flog, '''#############################################################
## Using TGAS catalog:'''
print>>flog, time.strftime('##The scripts runs at %Y-%m-%d %H:%M:%S\n',time.localtime(time.time()))
flog1= open('../logs/GDfit.tex', 'a')
print>>flog1, time.strftime('#The scripts runs at %Y-%m-%d %H:%M:%S',time.localtime(time.time()))
print>>flog1, '## result that can be used in latex.'

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
bar = 30
#bar = 1e10
print>>flog,  '## barrier: %d '%bar
print>>flog1, '## barrier: %d '%bar
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
## consider the systematic error in parallax of -0.25mas
#err = -0.25
#err = -0.25
err = 0.20
plx = plx*(1+err)
print>>flog1, ' ##plx = plx *( 1+ %.2f):'%err
print>>flog,  ' ##plx = plx *( 1+ %.2f):'%err
## Heliocentric distance. Unit: kpc.
r = 1.0/plx
## Galactic rectangle coordinate。 Unit：kpc.
X = r*cos(l)*cos(b)
Y = r*sin(l)*cos(b)
Z = r*sin(b)

## lower limit of r
rl = np.arange(0.0, 0.6, 0.1)
ru = 1.0
## limit of r: r >= 0.3kpc
for i in range(len(rl)):
    indr1 = np.where(r>=rl[i])
    indr2 = np.where(r<=ru)
    indr  = np.intersect1d(indr1, indr2)
    #print>>flog, '## limit of r: r > 0.3kpc:\n number   %d'%len(indr[0])
##     galactic height |z| <= 0.5kpc 
    indz = np.where(np.abs(Z)<=0.5) 
    #print>>flog, '## limit of z: |z| <= 0.5kpc:\n number   %d'%len(indz[0])
    ind = np.intersect1d(indr, indz)
    #print>>flog, '## Combine these two limitations:\n number   %d'%len(ind)
    
    ## extract the sub-sample.
    tags    = np.take(tag,    ind)
#    pmraTs  = np.take(pmraT,  ind)
#    pmdecTs = np.take(pmdecT, ind)
    pmraGs  = np.take(pmraG,  ind)
    pmdecGs = np.take(pmdecG, ind)
    raGs    = np.take(raG,    ind)
    decGs   = np.take(decG,   ind)
#    plxs    = np.take(plx,    ind)
    rs = np.take(r, ind)
    ls = np.take(l, ind)
    bs = np.take(b, ind)
#    Xs = np.take(X, ind)
#    Ys = np.take(Y, ind)           
#    Zs = np.take(Z, ind)
    
    ## proper motion from 
#    pmlTs, pmbTs, e_pmlTs, e_pmbTs = pmT(pmraTs, pmdecTs, np.ones(len(pmraTs)), np.ones(len(pmraTs)),\
#    raGs, decGs, ls)
    pmlGs, pmbGs, e_pmlGs, e_pmbGs = pmT(pmraGs, pmdecGs, np.ones(len(pmraGs)), np.ones(len(pmraGs)),\
    raGs, decGs, ls)
    
    print>>flog, "## %.1f kpc < r < %.1f kpc: TGAS"%(rl[i],ru)
    print>>flog1,"## %.1f kpc < r < %.1f kpc: TGAS"%(rl[i],ru)
    ParameterFit( pmlGs, pmbGs, ls, bs, rs, flog, flog1)
    ParameterFit1(pmlGs, pmbGs, ls, bs, rs, flog, flog1)
    ParameterFit2(pmlGs, pmbGs, ls, bs, rs, flog, flog1)
#    ParameterFit2(pmlGs, ls, bs, rs, flog)
    
#    print>>flog, "## %.1f kpc < r < 1.0kpc: TGAS-Tycho2"%rl[i]
#    print>>flog1,"## %.1f kpc < r < 1.0kpc: TGAS-Tycho2"%rl[i]
#    ParameterFit( pmlGs-pmlTs, pmbGs-pmbTs, ls, bs, rs, flog, flog1)
#    ParameterFit1(pmlGs-pmlTs, pmbGs-pmbTs, ls, bs, rs, flog, flog1)
#    ParameterFit2(pmlGs-pmlTs, ls, bs, rs, flog)
##print>>flog, '''\n################## consider the influence of systematic error 0.3mas
##'''
#### plx = plx + 0.3
##print>>flog, '1. plx = plx + 0.3:'
##ParameterFit(pmlGs, pmbGs, ls, bs, rPs, flog)
##ParameterFit1(pmlGs, pmbGs, ls, bs, rPs, flog)
##ParameterFit2(pmlGs, ls, bs, rPs, flog)
##print>>flog, '2. plx = plx - 0.3:'
#### plx = plx - 0.3
##ParameterFit(pmlGs, pmbGs, ls, bs, rMs, flog)
##ParameterFit1(pmlGs, pmbGs, ls, bs, rMs, flog)
##ParameterFit2(pmlGs, ls, bs, rMs, flog)
#
print'Done!'