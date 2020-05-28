# -*- coding: utf-8 -*-"""Created on Fri Mar 25 15:46:30 2016@author: NeoDivide the sources into 4 groups, each has the equal spherical area.And Normalized according to the apparent proper motionThis code is written for ranking sources by groupsOct 19, 2016: updated by Niu."""import numpy as npimport time## APM data file.apm = '../results/new_candidate4.apm'Sou = np.loadtxt(apm, usecols=(0,), dtype=str)V, V_RA, V_Dec, V_err, Va_E, Vd_E, RA, Dec\    = np.loadtxt(apm, usecols=range(1,9), unpack=True)dat1 = []dat2 = []dat3 = []dat4 = []for i in range(len(Sou)):    if Dec[i] < -30:        dat1.append((Sou[i], V[i], V_err[i],\        RA[i], V_RA[i], Va_E[i], Dec[i], V_Dec[i],Vd_E[i], V[i]/V_err[i]))    elif Dec[i] < 0:        dat2.append((Sou[i], V[i], V_err[i],\        RA[i], V_RA[i], Va_E[i], Dec[i], V_Dec[i],Vd_E[i], V[i]/V_err[i]))    elif Dec[i] < 30:        dat3.append((Sou[i], V[i], V_err[i],\        RA[i], V_RA[i], Va_E[i], Dec[i], V_Dec[i],Vd_E[i], V[i]/V_err[i]))    else :        dat4.append((Sou[i], V[i], V_err[i],\        RA[i], V_RA[i], Va_E[i], Dec[i], V_Dec[i],Vd_E[i], V[i]/V_err[i]))dtype = [('sou', 'S10'), ('u', float), ('u_e', float),\        ('alp', float), ('ua', float), ('ua_e', float),\        ('det', float), ('ud', float), ('ud_e', float), ('nu', float)]## sorted criterions = 'nu'#s = 'u'        a = np.array(dat1, dtype=dtype)        b1 = np.sort(a, order=[s]) a = np.array(dat2, dtype=dtype)        b2 = np.sort(a, order=[s]) a = np.array(dat3, dtype=dtype)        b3 = np.sort(a, order=[s]) a = np.array(dat4, dtype=dtype)        b4 = np.sort(a, order=[s]) [N1, N2, N3, N4] = [len(b1), len(b2), len(b3), len(b4)]## res_dir = '../results/'res = 'GR.rank'fou = open(res_dir + res, 'w')print>>fou, time.strftime('## This list was generated at %Y-%m-%d %H:%M:%S',\    time.localtime(time.time()))print>>fou, '''## original file: %s, sorted by: %s## data format: Source  pm  e_pm  RA  pmRA  e_pmRA  De  pmDe  e_pmDE##   deg for positions and uas/yr for proper motions##  source number of sub-groups: '''%(res, s), len(b1), len(b2), len(b3), len(b4) ##[N1, N2, N3, N4] = [77, 126, 161, 135] for 499sources#print len(b1), len(b2), len(b3), len(b4)i = 0while i < N3:    if i < N1:        b = np.sort([b1[i], b2[i], b3[i], b4[i]], order=[s])    elif i < N2:        b = np.sort([b2[i], b3[i], b4[i]], order=[s])    elif i < N4:        b = np.sort([b3[i], b4[i]], order=[s])    else:        b = [b3[i]]    for j in range(len(b)):        print >> fou, b[j][0] +"  %15.10f"*8\            %(b[j][1],b[j][2],b[j][3],b[j][4],b[j][5],b[j][6],b[j][7],b[j][8])        i += 1fou.close()print 'Done!'