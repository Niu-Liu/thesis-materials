# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 12:13:57 2016

@author: Neo

TRAnsformation:
Equatorial coordinate -->  Galactic coordinate

"""

import numpy as np
import sys
cos = np.cos
sin = np.sin
tan = np.tan
#cot = np.cot
pi  = np.pi
acos= np.arccos
asin= np.arcsin

## inclination
i = np.deg2rad(62.87175)     # RAd
## longtitude of ascending node Ω
l_an = np.deg2rad(32.93192)  # RAd
## the right ascension of N
RA_N = np.deg2rad(282.85948)

def cot(x):
    return cos(x)/sin(x)

def sec(x):
    return 1.0/cos(x)

def cord2arc(X):
    x, y, z = X/np.sqrt(np.sum(X**2))
#    if x**2+y**2+z**2 > 1:
#        print 'x^2+y^2+z^2>1. Illegal!'
##        sys.exit(1)
#    print x**2+y**2+z**2
        
    det = asin(z)
    cosa= x/cos(det)
    sina= y/cos(det)
    alpc= acos(cosa)
    if sina>0:
        alp = alpc
    else:
        alp = 2*pi - alpc
        
    return alp, det
    
def arc2cord(alp, det):
    '''
    alp/det: rad
    '''
    if not (0<alp<2*pi and -pi/2<det<pi/2):
        print 'Illegal! '
        sys.exit(1)
    
    x = cos(alp)*cos(det)
    y = sin(alp)*cos(det)
    z = sin(det)
    return [x, y, z]

# for test
# test1
#x = 0.2
#y = 0.3
#z = np.sqrt(1-x*x-y*y)
#print x, y, z
#a, b = cord2arc(x, y, z)
#print a, b
#x, y, z = arc2cord(a, b)
#print x, y, z
# test2
#b = pi/5
#a = pi*7/6
#print a, b
#x, y, z = arc2cord(a, b)     
#print x, y, z
#a, b = cord2arc(x, y, z)
#print a, b

def Ga2Eq_rect(XG):
    A = np.mat([[ -0.0548755604, -0.8734370902, -0.4838350155], \
                [ +0.4941094279, -0.4448296300, +0.7469822445], \
                [ -0.8676661490, -0.1980763734, +0.4559837762]])
    XE = np.dot(np.transpose(A), XG)
    
    return XE
    
def Eq2Ga_rect(XE):
    A = np.mat([[ -0.0548755604, -0.8734370902, -0.4838350155], \
                [ +0.4941094279, -0.4448296300, +0.7469822445], \
                [ -0.8676661490, -0.1980763734, +0.4559837762]])
    XG = np.dot(A, XE)
    
    return XG

def Eq2Ga(RA, DE):
    '''
    RA/DE: rad --> l/b: rad
    S(l,b)=A*S(RA,DE)
    A = [ [ -0.0548755604, -0.8734370902, -0.4838350155]
          [ +0.4941094279, -0.4448296300, +0.7469822445]
          [ -0.8676661490, -0.1980763734, +0.4559837762]]
    '''
    N = len(RA)
    
    A = np.mat([[ -0.0548755604, -0.8734370902, -0.4838350155], \
          [ +0.4941094279, -0.4448296300, +0.7469822445], \
          [ -0.8676661490, -0.1980763734, +0.4559837762]])
    l = []
    b = []
          
    for i in range(N):
        XE = np.array(arc2cord(RA[i], DE[i]))
        XG = np.dot(A, XE)
        l0, b0 = cord2arc(np.array(XG)[0])
        l.append(l0)
        b.append(b0)
        
    return np.array(l), np.array(b)
        

def CalPhi(RA, DE, l):
    '''
    \phi is the included angle. And the value is given by the following equations.
    cot(\phi) = cos(i)*cos(DE)*sec(RA-RA_N) + sin(DE)*tan(RA-RA_N)
    sin(\phi) = sin(i)*cos(l-l_an)/cos(DE)
    
    NB: units are rads for RA, DE and l
    '''
    cotPhi = cos(i)*cos(DE)*sec(RA-RA_N) + sin(DE)*tan(RA-RA_N)
    sinPhi = sin(i)*cos(l-l_an)/cos(DE)
    cosPhi = cotPhi*sinPhi
#    cosPhi = 
    
    return sinPhi, cosPhi
    
def pmT(pmRA, pmDE, e_pmRA, e_pmDE, RA, DE, l):
    '''
    pmLon = +pmRA*cos(\phi) + pmDE*sin(\phi)
    pmLat = -pmRA*sin(\phi) + pmDE*cos(\phi)
    
    NB: pmRA = mu_RA*cos(DE), pmLon = mu_lon*cos(b)
    '''
    sinPhi, cosPhi = CalPhi(RA, DE, l)
    pmLon = +pmRA*cosPhi + pmDE*sinPhi
    pmLat = -pmRA*sinPhi + pmDE*cosPhi
    
    e_pmLon = np.sqrt(e_pmRA**2*cosPhi**2 + e_pmDE**2*sinPhi**2)
    e_pmLat = np.sqrt(e_pmRA**2*sinPhi**2 + e_pmDE**2*cosPhi**2)
    
    return pmLon, pmLat, e_pmLon, e_pmLat