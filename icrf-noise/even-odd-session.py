#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: even-odd-session.py
"""
Created on Wed May 30 15:13:12 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Sep 29 2018: This code will not be used any longer. The replacement is decimation_solution_difference.ipnb.

"""

import numpy as np
from numpy import cos, deg2rad, sqrt
import sys
import time
# My modules
# from rewrite_sou import rewrite_sou
from nor_sep import pos_max_calc, overall_err_calc
from list_crossmatch import list_crossmatch
from error_plot import (
    error_vs_numses, error_vs_numses2,
    error_vs_numobs, error_vs_numobs2, maxerror_vs_numses,
    maxerror_vs_numobs, maxerror_vs_numses2, maxerror_vs_numobs2,
    overallerror_vs_numses, overallerror_vs_numobs,
    overallerror_vs_numses2, overallerror_vs_numobs2,
    comparison_plot, comparison_plot2, comparison_plot3, comparison_plot3_1,
    difference_vs_error_plot, sf_nf_dec_plot)

# -----------------------------  FUNCTIONS -----------------------------

# --------------------------------- MAIN -------------------------------
# Rewrite catalogs
# print("Rewrite .sou file:")
# # the solution of odd sessions
# rewrite_sou("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
#             "opa2018b-ga15-odd/opa2018b-ga15-odd.sou")
# # the solution of even sessions
# rewrite_sou("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
#             "opa2018b-ga15-even/opa2018b-ga15-even.sou")

lab = sys.argv[3]

# Load data
# odd sessions
print("Odd session:")
# cat1 = "../data/opa2018b-ga15-odd.cat"
cat1 = "../data/%s" % sys.argv[1]
sou1 = np.genfromtxt(cat1, dtype=str, usecols=(0,))
RA1, Dec1, RAc_err1, Dec_err1, corr1 = np.genfromtxt(
    cat1, usecols=range(2, 7), unpack=True)
num_ses1, num_obs1 = np.genfromtxt(
    cat1, usecols=range(10, 12), dtype=int, unpack=True)

print("Median formal error (uas):\n"
      "  RA:  %6.1f\n"
      "  Dec: %6.1f" % (np.median(RAc_err1) * 1.e3,
                        np.median(Dec_err1) * 1.e3))
print("range of N_session: [%d, %d]" % (min(num_ses1), max(num_ses1)))
print("range of N_observation: [%d, %d]" % (min(num_obs1), max(num_obs1)))

# ellipe semi-major axis
sig_pos_max1 = pos_max_calc(RAc_err1, Dec_err1, corr1)

# overall formal uncertainty
overall_err1 = overall_err_calc(RAc_err1, Dec_err1, corr1)

# Error plot for odd solution
error_vs_numses(RAc_err1, Dec_err1, num_ses1,
                "../plots/odd-err-Nses-%s.eps" % lab)
error_vs_numobs(RAc_err1, Dec_err1, num_obs1,
                "../plots/odd-err-Nobs-%s.eps" % lab)

maxerror_vs_numses(sig_pos_max1, num_ses1,
                   "../plots/odd-maxerr-Nses-%s.eps" % lab)
maxerror_vs_numobs(sig_pos_max1, num_obs1,
                   "../plots/odd-maxerr-Nobs-%s.eps" % lab)

overallerror_vs_numses(overall_err1, num_ses1,
                       "../plots/odd-overallerr-Nses-%s.eps" % lab)
overallerror_vs_numobs(overall_err1, num_obs1,
                       "../plots/odd-overallerr-Nobs-%s.eps" % lab)


# even sessions
print("Even session:")
cat2 = "../data/%s" % sys.argv[2]
sou2 = np.genfromtxt(cat2, dtype=str, usecols=(0,))
RA2, Dec2, RAc_err2, Dec_err2, corr2 = np.genfromtxt(
    cat2, usecols=range(2, 7), unpack=True)
num_ses2, num_obs2 = np.genfromtxt(cat2, usecols=range(10, 12),
                                   dtype=int, unpack=True)

print("Median formal error (uas):\n"
      "  RA:  %6.1f\n"
      "  Dec: %6.1f" % (np.median(RAc_err2) * 1.e3,
                        np.median(Dec_err2) * 1.e3))
print("range of N_session: [%d, %d]" % (min(num_ses2), max(num_ses2)))
print("range of N_observation: [%d, %d]" % (min(num_obs2), max(num_obs2)))

# ellipe semi-major axis
sig_pos_max2 = pos_max_calc(RAc_err2, Dec_err2, corr2)

# overall formal uncertainty
overall_err2 = overall_err_calc(RAc_err2, Dec_err2, corr2)

# Error plot for even solution
error_vs_numses(RAc_err2, Dec_err2, num_ses2,
                "../plots/even-err-Nses-%s.eps" % lab)
error_vs_numobs(RAc_err2, Dec_err2, num_obs2,
                "../plots/even-err-Nobs-%s.eps" % lab)

maxerror_vs_numses(sig_pos_max2, num_ses2,
                   "../plots/even-maxerr-Nobs-%s.eps" % lab)
maxerror_vs_numobs(sig_pos_max2, num_obs2,
                   "../plots/even-maxerr-Nses-%s.eps" % lab)

overallerror_vs_numses(overall_err2, num_ses2,
                       "../plots/even-overallerr-Nses-%s.eps" % lab)
overallerror_vs_numobs(overall_err2, num_obs2,
                       "../plots/even-overallerr-Nobs-%s.eps" % lab)


# Plot for both solution
error_vs_numses2(RAc_err1, Dec_err1, num_ses1, RAc_err2, Dec_err2, num_ses2,
                 "../plots/odd-even-err-Nses-%s.eps" % lab)
error_vs_numobs2(RAc_err1, Dec_err1, num_obs1, RAc_err2, Dec_err2, num_obs2,
                 "../plots/odd-even-err-Nobs-%s.eps" % lab)

maxerror_vs_numses2(sig_pos_max1, num_ses1, sig_pos_max2, num_ses2,
                    "../plots/odd-even-maxerr-Nses-%s.eps" % lab)
maxerror_vs_numobs2(sig_pos_max1, num_obs1, sig_pos_max2, num_obs2,
                    "../plots/odd-even-maxerr-Nobs-%s.eps" % lab)

overallerror_vs_numses2(overall_err1, num_ses1, overall_err2, num_ses2,
                        "../plots/odd-even-overallerr-Nses-%s.eps" % lab)
overallerror_vs_numobs2(overall_err1, num_obs1, overall_err2, num_obs2,
                        "../plots/odd-even-overallerr-Nobs-%s.eps" % lab)

# Cross-match between between two solutions
comsou, ind1, ind2 = list_crossmatch(sou1, sou2)

# Verify the result of cross-match
soucom1 = np.take(sou1, ind1)
soucom2 = np.take(sou2, ind2)
for i, (comsoui, soucom1i, soucom2i) in enumerate(
        zip(comsou, soucom1, soucom2)):

    if comsoui != soucom1i:
        print("%dth source %s are not consistented in list1 %s." %
              (i, comsoui, soucom1i))

    if comsoui != soucom2i:
        print("%dth source %s are not consistented in list2 %s." %
              (i, comsoui, soucom2i))


# Extract data
# odd sessions
RA_com1 = np.take(RA1, ind1)
Dec_com1 = np.take(Dec1, ind1)
RAc_err_com1 = np.take(RAc_err1, ind1)
Dec_err_com1 = np.take(Dec_err1, ind1)
num_ses_com1 = np.take(num_ses1, ind1)
num_obs_com1 = np.take(num_obs1, ind1)
sig_pos_max_com1 = np.take(sig_pos_max1, ind1)

# even sessions
RA_com2 = np.take(RA2, ind2)
Dec_com2 = np.take(Dec2, ind2)
RAc_err_com2 = np.take(RAc_err2, ind2)
Dec_err_com2 = np.take(Dec_err2, ind2)
num_ses_com2 = np.take(num_ses2, ind2)
num_obs_com2 = np.take(num_obs2, ind2)
sig_pos_max_com2 = np.take(sig_pos_max2, ind2)


# # Plot for error comparsion
comparison_plot(RAc_err_com1, RAc_err_com2,
                "$\\sigma_{\\alpha^*}$ (mas) - odd",
                "$\\sigma_{\\alpha^*}$ (mas) - even",
                "../plots/odd-even-com-RAcerr-%s.eps" % lab)
comparison_plot(Dec_err_com1, Dec_err_com2,
                "$\\sigma_\\delta$ (mas) - odd",
                "$\\sigma_\\delta$ (mas) - even",
                "../plots/odd-even-com-Decerr-%s.eps" % lab)
comparison_plot(num_ses_com1, num_ses_com2,
                "$N_{session}$ - odd",
                "$N_{session}$ - even",
                "../plots/odd-even-com-Numses-%s.eps" % lab)
comparison_plot(num_obs_com1, num_obs_com2,
                "$N_{observation}$ - odd",
                "$N_{observation}$ - even",
                "../plots/odd-even-com-Numobs-%s.eps" % lab)
comparison_plot(sig_pos_max_com1, sig_pos_max_com2,
                "Elliptical error (mas) - odd",
                "Elliptical error (mas) - even",
                "../plots/odd-even-com-maxerr-%s.eps" % lab)


# Calculate the differences (even - odd)
# deg -> mas
dRAc = (RA_com2 - RA_com1) * cos(deg2rad(Dec_com2)) * 3.6e6
dDec = (Dec_com2 - Dec_com1) * 3.6e6
dRAc_err_2 = RAc_err_com1**2 + RAc_err_com2**2
dRAc_err = sqrt(dRAc_err_2)
dDec_err_2 = Dec_err_com1**2 + Dec_err_com2**2
dDec_err = sqrt(dDec_err_2)
num_ses_mean = (num_ses_com1 + num_ses_com2) / 2.
num_obs_mean = (num_obs_com1 + num_obs_com2) / 2.


# print the result.
fout = open("../data/odd-even-%s.cat_diff" % lab, "w")

# Header
print("# VLBI Celestial Reference Frame Solution OPA2018b decimation test\n"
      "#\n"
      "# Columns  Units   Meaning\n"
      "#    1     --      IVS designation\n"
      "#    2     deg     Right ascension(odd session)\n"
      "#    3     deg     Declination(odd session)\n"
      "#    4     mas     Difference of RA*cos(Dec)\n"
      "#    5     mas     Difference of Declination\n"
      "#    6     mas     Formal uncertainty of difference in "
      "RA*cos(Dec)\n"
      "#    7     mas     Formal uncertainty of difference in "
      "declination\n"
      "#    8     --      Mean number of sessions in two solutions\n"
      "#    9     --      Mean number of observations in two solutions\n"
      "# Created date: %s\n#"
      % time.strftime("%d/%m/%Y", time.localtime()), file=fout)

for (comsoui, RA_com1i, Dec_com1i, dRAci, dDeci, dRAc_erri, dDec_erri,
     num_ses_meani, num_obs_meani) in zip(
        comsou, RA_com1, Dec_com1, dRAc, dDec, dRAc_err, dDec_err,
        num_ses_mean, num_obs_mean):
    print("%-8s  %+10.4f  %+10.4f  %+8.4f  %+8.4f  %8.4f  %8.4f  %6d  %8d" %
          (comsoui, RA_com1i, Dec_com1i, dRAci, dDeci, dRAc_erri, dDec_erri,
           num_ses_meani, num_obs_meani), file=fout)
fout.close()


# Plot for differences vs error
difference_vs_error_plot(dRAc_err_2, dRAc**2,
                         "$\\sigma_{\\Delta\\alpha^*}^2 (mas^2)$",
                         "$\\Delta^2_{\\alpha^*} (mas^2)$",
                         "../plots/odd-even-dRAc-dRAcerr-%s.eps" % lab)
difference_vs_error_plot(dDec_err_2, dDec**2,
                         "$\\sigma_{\\Delta\\delta}^2 (mas^2)$",
                         "$\\Delta^2_{\\delta} (mas^2)$",
                         "../plots/odd-even-dDec-dDecerr-%s.eps" % lab)
difference_vs_error_plot(dRAc_err_2, dRAc,
                         "$\\sigma_{\\Delta\\alpha^*}^2 (mas^2)$",
                         "$\\Delta^2_{\\alpha^*} (mas^2)$",
                         "../plots/odd-even-dRAc-dRAcerr-%s.eps" % lab)
difference_vs_error_plot(dDec_err_2, dDec**2,
                         "$\\sigma_{\\Delta\\delta}^2 (mas^2)$",
                         "$\\Delta^2_{\\delta} (mas^2)$",
                         "../plots/odd-even-dDec-dDecerr-%s.eps" % lab)
# --------------------------------- END --------------------------------
