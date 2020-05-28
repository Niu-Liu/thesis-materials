# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 10:36:34 2016

@author: Neo

cp the figure ../plot/ into ../NLiu2016/
"""

import os

## figure names in ../plot/
l1 = ['Observation_span.eps',    \
      'Number_of_session.eps',   \
      'Observation_history.eps', \
      'DEV_plot.eps',            \
      'LinearDriftOf4Sets.eps',  \
      'Rotation_No.eps',         \
      'RotAndGli.eps',           \
      'Simulation.eps',          \
      'Quality.eps',             \
      'LD_ours.eps',             \
      '39Special.eps',           \
      '294SouDis.eps',           \
      '323SouDis.eps'            \
      ]
      
#       WADEV_de.eps                  ld_4.eps
#331SouDis.eps                 Linear_drift_4.eps            Observation_span.eps          WADEV_ra.eps                  rot_num.eps
#Combination.eps               Linear_drift_icrf1.eps        Quality.eps                   WDEV_de.eps                   seslen.eps
#                  Linear_drift_icrf2.eps        RotAndGli.eps                 WDEV_ra.eps
#GRank_rot.eps                 Linear_drift_mfv247.eps       Rot_Gli.eps                   gli_num.eps
#GrRank_rot.eps                Lineardrift2error_special.eps Rotation_No.eps               group_num.eps
#LD_ours.eps                            Simulation.eps                grp_gli_num.eps']

l2 = ['fig1\(a\).eps', \
      'fig1\(b\).eps', \
      'fig2.eps',    \
      'fig3.eps',    \
      'fig4.eps',    \
      'fig5.eps',    \
      'fig6.eps',    \
      'fig7.eps',    \
      'fig8.eps',    \
      'fig9.eps',    \
      'fig11.eps',   \
      'fig10\(a\).eps',\
      'fig10\(b\).eps',\
      ]
      
if len(l1) != len(l2):
    print('Error! Unequal length of these two lists.')
else:
    for i in range(len(l1)):
        os.system('cp ../plot/'+l1[i]+' ../manuscript/figures/'+l2[i])
        print('copy figure ' + l1[i] + ' : Done!')
        
print('All Done!')