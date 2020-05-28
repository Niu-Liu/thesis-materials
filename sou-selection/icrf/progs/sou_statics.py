# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 22:10:03 2015

@author: Neo
"""

# time series statics of sources

# initial value of variables
cat='../catalog/OPA.cat'  #catalog file, read
sta='../results/sou.sta'  #statics file, write

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
            
# 
#       
fsta=open(sta,'w')
fsta.write('#---------------------------------------------------------\n'+\
'#sou_name    ses_num    obs_begin    obs_end    obs_span(yr)\n'+\
'#--------------------------------------------------------------------\n')
fcat=open(cat,'r')
con=fcat.readline()
while len(con):
    souf=con.strip('\n')
    fdat=open('../data/opa/'+souf +'.dat','r')
    datcon=fdat.readlines()
    fdat.close()

#source data format
#Obs_moment   R.A.(o)   Dec.(o)  Err. (mas)  Corr. Del.  IERS   IVS
    seslen=len(datcon)-3
#    soun=souf[0:8]
    dat_b=str(datcon[3]).split()
    dat_e=str(datcon[-1]).split()
    obs_b=ADepo(float(dat_b[0]))
    obs_e=ADepo(float(dat_e[0]))
    fsta.write(souf+'    '+str(seslen)+'    '+str(obs_b)+'    '+\
    str(obs_e)+'    '+str(obs_e-obs_b)+'\n')
    
    con=fcat.readline()
    
fcat.close()

#code written before
#for i in range(len(catcon)):
#    souf=str(catcon[i]).strip('\n')
#    fdat=open(souf,'r')
#    datcon=fdat.readlines()
#    fdat.close()
#    #source data format
#    #Obs_moment   R.A.(o)   Dec.(o)  Err. (mas)  Corr. Del.  IERS   IVS
#    seslen=len(datcon)-3
#    soun=souf[0:8]
#    dat_b=str(datcon[3]).split()
#    dat_e=str(datcon[-1]).split()
#    obs_b=ADepo(float(dat_b[0]))
#    obs_e=ADepo(float(dat_e[0]))
#    fsta.write(soun+'    '+str(seslen)+'    '+str(obs_b)+'    '+\
#    str(obs_e)+'    '+str(obs_e-obs_b)+'\n')

fsta.close()
print 'Done!'    
