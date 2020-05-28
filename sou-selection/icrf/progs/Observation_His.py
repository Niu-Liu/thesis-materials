# -*- coding: utf-8 -*-
"""
Created on Thu  Mar 3 16:36:40 2016

@author: Neo
"""

# For session coordinate time series, plot the session_time vs declination
#

##Some Function to be used later.
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
        
#convert Julian epoch to A.D. epoch
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
        elif date >= day and date <= con_list[i+1][1]:
            return yr+(date-day)/leap_year(int(yr))
            
import matplotlib.pyplot as plt
import numpy as np

#ICRF2 defining sources data file
#icrf2 = '../catalog/icrf2-defining.dat'
#defn = []
#deff = open(icrf2)
#conf = deff.readlines()
#for i in range(20,len(conf)):
#    dlin=str(conf[i]).split()
#    defn.append(dlin[2])
#deff.close()

## For candidates
#cat='../catalog/new_candidate2.cat'  #catalog file, read
##
cat='../catalog/ex_def.cat'
fcat = open(cat,'r')
con=fcat.readline()

plt.figure(figsize=(20,15))

while len(con):
    souf=con.strip('\n')
    fdat=open('../data/opa/'+souf +'.dat','r')

#data format
#Obs_moment   R.A.(o)   Dec.(o)  Err. (mas)  Corr. Del.  IERS   IVS
    t_ses = []
    Dec_ses = []
    datcon = fdat.readlines()
    datcon = datcon[3:]
    for i in range(len(datcon)):
        dat = str(datcon[i]).split()
        t_ses.append(ADepo(float(dat[0])))
        Dec_ses.append(float(dat[2]))
        
    Dec = sum(Dec_ses) / len(Dec_ses)
    dec = np.ones(len(t_ses)) * Dec
    
    print len(t_ses)
    
#    if souf in defn:
#        plt.plot(t_ses, dec, 'r.', markersize=5)
#    else:
#        plt.plot(t_ses, dec, 'b.', markersize=5)
        
    plt.plot(t_ses, dec, 'b.', markersize=5)
                    
    fdat.close()
    con=fcat.readline()

    
plt.xlim([1979.0, 2016.0])   
plt.ylim([-90.0, 90.0])     
plt.yticks([-90, -60, -30, 0, 30, 60 ,90], \
    ['-90', '-60', '-30', '0', '30', '60', '90'], fontsize=30)
plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015], \
    ['1980', '', '1990', '', '2000', '', '2010', ''], fontsize=30)    
plt.ylabel('Declination($\degree$)', fontsize=35)
#plt.savefig('../plot/Observation_history.png')
#plt.savefig('../plot/Observation_history.pdf')
plt.show()
#plt.close()

#fcat.close()
#foutp.close()
print 'Done!'    
