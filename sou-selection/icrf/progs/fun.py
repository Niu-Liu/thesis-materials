# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 14:01:32 2015

@author: Neo
"""
#!/bin/python
#function defined by myself.

import numpy as np
from math import *
#pi=math.pi
#cos=math.cos
#sin=math.sin
#sqrt=math.sqrt

#Leap year
def leap_year(year):
    flag=0
    if year%100:
        if year%4:
            flag=1
    elif year%400 :
        flag=1
        
    if flag:
        return 365
    else:
        return 366
        
#convert MJD(Modified Julian Date) epoch to A.D. epoch
def ADepo(date):
    con_list=[\
    [1979.0,43874.0],\
    [1980.0,44240.0],\
    [1981.0,44605.0],\
    [1982.0,44970.0],\
    [1983.0,45335.0],\
    [1984.0,45700.0],\
    [1985.0,46066.0],\
    [1986.0,46431.0],\
    [1987.0,46796.0],\
    [1988.0,47161.0],\
    [1989.0,47527.0],\
    [1990.0,47892.0],\
    [1991.0,48257.0],\
    [1992.0,48622.0],\
    [1993.0,48988.0],\
    [1994.0,49353.0],\
    [1995.0,49718.0],\
    [1996.0,50083.0],\
    [1997.0,50449.0],\
    [1998.0,50814.0],\
    [1999.0,51179.0],\
    [2000.0,51544.0],\
    [2001.0,51910.0],\
    [2002.0,52275.0],\
    [2003.0,52640.0],\
    [2004.0,53005.0],\
    [2005.0,53371.0],\
    [2006.0,53736.0],\
    [2007.0,54101.0],\
    [2008.0,54466.0],\
    [2009.0,54832.0],\
    [2010.0,55197.0],\
    [2011.0,55562.0],\
    [2012.0,55927.0],\
    [2013.0,56293.0],\
    [2014.0,56658.0],\
    [2015.0,57023.0]\
    ]
    for i in range(len(con_list)):
        yr=con_list[i][0]
        day=con_list[i][1]
        if i==len(con_list)-1:
            return yr+(date-day)/leap_year(int(yr))
        elif date>=day and date<=con_list[i+1][1]:
            return yr+(date-day)/leap_year(int(yr))
            
##Apparent proper motion
#Fitting curves
#alp = mu_a*(t-t0)+alp0
#det = mu_d*(t-t0)+det0

def linlsq(alp,siga,det,sigd,epo):
    L=len(epo)
    F=np.empty([L,2], dtype=float)
    ya=np.empty([L,1], dtype=float)
    yd=np.empty([L,1], dtype=float)
    Pa=np.zeros([L,L], dtype=float)
    Pd=np.zeros([L,L], dtype=float)
    
    for i in range(len(epo)):
        F[i,0]=epo[i]
        F[i,1]=1
        Pa[i,i]=1.0/siga[i]**2
        Pd[i,i]=1.0/sigd[i]**2
        ya[i,0]=alp[i]
        yd[i,0]=det[i]
    
    F=np.matrix(F)
    ya=np.matrix(ya)
    yd=np.matrix(yd)
    Pa=np.matrix(Pa)
    Pd=np.matrix(Pd)
    
    Aa=(F.T.dot(Pa)).dot(F)
    Ad=(F.T.dot(Pd)).dot(F)
    
    Va=Aa.I
    Vd=Ad.I
    
    ca=Va.dot(F.T).dot(Pa).dot(ya)
    cd=Vd.dot(F.T).dot(Pd).dot(yd)
    
    ua=ca[0,0]
    alp0=ca[1,0]/3.6e9
    ud=cd[0,0]
    det0=cd[1,0]/3.6e9
    
    sigua=sqrt(Va[0,0])
    sigud=sqrt(Vd[0,0])
    siga0=sqrt(Va[1,1]/3.6e9)
    sigd0=sqrt(Vd[1,1]/3.6e9)
    
    chia=(ya-F.dot(ca)).T.dot(Pa).dot(ya-F.dot(ca))[0,0]
    chid=(yd-F.dot(cd)).T.dot(Pd).dot(yd-F.dot(cd))[0,0]
    
    return [ua,ud,sigua,sigud,alp0,det0,siga0,sigd0,chia,chid]

def pm(alp,siga,det,sigd,epo,tb):
    i=0
    
    #tb stands for the beginning epoch of +++++++++++`````````````````````````````````````````````
    while i<len(epo):
        if epo[i]<tb or siga==0.0 or sigd==0.0:
            del epo[i],alp[i],siga[i],det[i],sigd[i]
        else:
            i+=1

    [ua,ud,sigua,sigud,alp0,det0,siga0,sigd0,chia,chid]\
    =linlsq(alp,siga,det,sigd,epo)
    
    u=math.sqrt((ua*cos(det0/180*pi))**2+ud**2)
    sigu=math.sqrt((sigua*cos(det0/180*pi))**2+sigud**2)
    
    return [u,sigu,ua,ud,sigua,sigud,alp0,det0,siga0,sigd0,chia,chid]
    
def axis_rot(ua,ud,sigua,sigud,alp0,det0):
    L=len(ua)
    Fa=np.empty([L,3], dtype=float)
    Fd=np.empty([L,3], dtype=float)
    ya=np.empty([L,1], dtype=float)
    yd=np.empty([L,1], dtype=float)
    Pa=np.zeros([L,L], dtype=float)
    Pd=np.zeros([L,L], dtype=float)
    
    for i in range(L):
        alp=radians(alp0[i])
        det=radians(det0[i])
        
        Fa[i,0] = cos(alp)*sin(det)
        Fa[i,1] = sin(alp)*sin(det)
        Fa[i,2] =-cos(det)
        
        Fd[i,0] =-sin(alp)
        Fd[i,1] = cos(alp)
        Fd[i,2] = 0.0      
        
        Pa[i,i] = 1.0/sigua[i]**2
        Pd[i,i] = 1.0/sigud[i]**2
        
        ya[i,0] = ua[i]
        yd[i,0] = ud[i]
    
    Fa=np.matrix(Fa)
    Fd=np.matrix(Fd)
    ya=np.matrix(ya)
    yd=np.matrix(yd)
    Pa=np.matrix(Pa)
    Pd=np.matrix(Pd)
    
    Aa=(Fa.T.dot(Pa)).dot(Fa)
    Ad=(Fd.T.dot(Pd)).dot(Fd)
    A=Aa+Ad
    
    V=A.I
    
    c=V.dot((Fa.T).dot(Pa).dot(ya)+(Fd.T).dot(Pd).dot(yd))
    
    wx=float(c[0])
    wy=float(c[1])
    wz=float(c[2])
    
    w=sqrt(wx**2+wy**2+wz**2)
    
    sigwx=sqrt(V[0,0])
    sigwy=sqrt(V[1,1])
    sigwz=sqrt(V[2,2])
    
    sigw=sqrt(sigwx**2+sigwy**2+sigwz**2)
    
    sigwxy=V[0,1]/sigwx/sigwy
    sigwxz=V[0,2]/sigwx/sigwz
    sigwyz=V[1,2]/sigwy/sigwz

    chia=(ya-Fa.dot(c)).T.dot(Pa).dot(ya-Fa.dot(c))[0,0]
    chid=(yd-Fd.dot(c)).T.dot(Pd).dot(yd-Fd.dot(c))[0,0]
    chi=chia+chid
    
    return [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,sigwxy,sigwxz,sigwyz,chi]
    
def Rotation_Glide(ua,ud,sigua,sigud,alp0,det0):
    L=len(ua)
    Fa=np.empty([L,6], dtype=float)
    Fd=np.empty([L,6], dtype=float)
    ya=np.empty([L,1], dtype=float)
    yd=np.empty([L,1], dtype=float)
    Pa=np.zeros([L,L], dtype=float)
    Pd=np.zeros([L,L], dtype=float)
    
    for i in range(L):
        alp=radians(alp0[i])
        det=radians(det0[i])
        
        Fa[i,0] = cos(alp)*sin(det)
        Fa[i,1] = sin(alp)*sin(det)
        Fa[i,2] =-cos(det)
        Fd[i,3] =-sin(alp)
        Fd[i,4] = cos(alp)
        Fd[i,5] = 0.0
        
        Fd[i,0] =-sin(alp)
        Fd[i,1] = cos(alp)
        Fd[i,2] = 0.0
        Fd[i,3] =-cos(alp)*sin(det)
        Fd[i,4] =-sin(alp)*sin(det)
        Fd[i,5] = cos(det)
        
        Pa[i,i] = sigua[i]**-2
        Pd[i,i] = sigud[i]**-2
        
        ya[i,0]=ua[i]
        yd[i,0]=ud[i]
    
    Fa=np.matrix(Fa)
    Fd=np.matrix(Fd)
    ya=np.matrix(ya)
    yd=np.matrix(yd)
    Pa=np.matrix(Pa)
    Pd=np.matrix(Pd)
    
    Aa=(Fa.T.dot(Pa)).dot(Fa)
    Ad=(Fd.T.dot(Pd)).dot(Fd)
    A=Aa+Ad
    
    V=A.I
    y=(Fa.T).dot(Pa).dot(ya)+(Fd.T).dot(Pd).dot(yd)    
    
    c=V.dot(y)
    
    wx=float(c[0])
    wy=float(c[1])
    wz=float(c[2])
    
    gx=float(c[3])
    gy=float(c[4])
    gz=float(c[5])
    
    w=sqrt(wx**2+wy**2+wz**2)
    g=sqrt(gx**2+gy**2+gz**2)
    
    sigwx = sqrt(V[0,0])
    sigwy = sqrt(V[1,1])
    sigwz = sqrt(V[2,2])
    
    siggx = sqrt(V[3,3])
    siggy = sqrt(V[4,4])
    siggz = sqrt(V[5,5])
    
    sigw = sqrt(sigwx**2+sigwy**2+sigwz**2)
    sigg = sqrt(siggx**2+siggy**2+siggz**2)

#    sig = [sigwx, sigwy, sigwz, siggx, siggy, siggz ]
#    
#    V0 = V
#    for i in range(6):
#        j = i
#        while j < 6:
#            V0[i,j] = V0[i,j]/sig[i]/sig[j]
#            j += 1
    
#    print V0        
#    

    de0 = atan(gz/sqrt(gx**2+gy**2))
    ra0 = atan(gy/gx) + pi
#    if gy/gx > 0:
#        if gx/cos(de0) <0:
#            ra0 += pi
#    else:
#        if gx/cos(de0) <0:
#            ra0 += pi
#        else:
#            ra0 += 2*pi
    print degrees(ra0), degrees(de0)

    chia=(ya-Fa.dot(c)).T.dot(Pa).dot(ya-Fa.dot(c))[0,0]
    chid=(yd-Fd.dot(c)).T.dot(Pd).dot(yd-Fd.dot(c))[0,0]
    chi=chia+chid
    
    return [w,sigw,wx,wy,wz,sigwx,sigwy,sigwz,\
            g,sigg,gx,gy,gz,siggx,siggy,siggz,chi]
    
def sou_division(det,N):
    '''
    return the (N-2) nodes  when dividing sources into N groups.
    '''
    num=len(det)/N
    dec_div=np.empty(N-2)
    
    det_sort=np.sort(det)
    for i in range(N-2):
        dec_div[i]=det_sort[num*(i+1)-1]
        
    return dec_div
    
def mean_apm(sou,defsou,u,sigu):
    f=0.0
    p=0.0
    for i in range(len(sou)):
        if sou[i] in defsou:
            pi=1.0/sigu[i]**2
            fi=u[i]*pi
            
            f+=fi
            p+=pi
            
    return f/p

def mean_dec(sou,det0):
#   Calculate the mean declination of source ensemble.
    N=len(sou)
    D=0.0    
    for i in range(N):
        D+=det0[i]
        
    return D/N
    
def simulation(alp, det):
#    For GC
#    alp0 = radians(266.25)
#    det0 = radians(-29)
    alp0 = radians(0)
    det0 = radians(90)
    A = 5
    
    d1 = A*cos(alp0)*cos(det0)
    d2 = A*sin(alp0)*cos(det0)
    d3 = A*sin(det0)
    
    ua = [-d1*sin(radians(alp[i])) + d2*cos(radians(alp[i])) \
         for i in range(len(alp))]
    ud = [-d1*cos(radians(alp[i]))*sin(radians(det[i]))-\
    d2*sin(radians(alp[i]))*sin(radians(det[i]))+d3*cos(radians(det[i]))\
         for i in range(len(alp))]
             
    return [ua, ud]
