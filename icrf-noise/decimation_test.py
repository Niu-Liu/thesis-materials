#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: VLBI_internal_error.py
"""
Created on Thu May 17 11:31:14 2018

@author: Neo(liuniu@smail.nju.edu.cn)


This code is used to determine the scaling factor and noise floor of the
VLBI formal error.

Sep 29 2018: This code will not be used any longer. The replacement is decimation_solution_difference.ipnb.

"""

import numpy as np
from numpy import cos, deg2rad, sqrt
from scipy.optimize import curve_fit
import time
# My modules
from nor_sep import pos_max_calc, overall_err_calc
from list_crossmatch import list_crossmatch
from error_plot import comparison_plot, difference_vs_error_plot


# -----------------------------  FUNCTIONS -----------------------------
# --------------------------------- MAIN -------------------------------
# Load data
# odd sessions
# dat_dir1 = ("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
#             "opa2018b-ga15-odd")

cat1 = "%s/opa2018b-ga15-odd.cat" % dat_dir1
sou1 = np.genfromtxt(cat1, dtype=str, usecols=(0,))
RA1, Dec1, RAc_err1, Dec_err1, corr1 = np.genfromtxt(
    cat1, usecols=range(2, 7), unpack=True)
num_ses1, num_obs1 = np.genfromtxt(cat1, usecols=range(10, 12),
                                   dtype=int, unpack=True)


# ellipe semi-major axis
sig_pos_max1 = pos_max_calc(RAc_err1, Dec_err1, corr1)

# overall formal uncertainty
overall_err1 = overall_err_calc(RAc_err1, Dec_err1, corr1)


# even sessions
dat_dir2 = ("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
            "opa2018b-ga15-even")

cat2 = "%s/opa2018b-ga15-even.cat" % dat_dir2
sou2 = np.genfromtxt(cat2, dtype=str, usecols=(0,))
RA2, Dec2, RAc_err2, Dec_err2, corr2 = np.genfromtxt(
    cat2, usecols=range(2, 7), unpack=True)
num_ses2, num_obs2 = np.genfromtxt(cat2, usecols=range(10, 12),
                                   dtype=int, unpack=True)


# ellipe semi-major axis
sig_pos_max2 = pos_max_calc(RAc_err2, Dec_err2, corr2)

# overall formal uncertainty
overall_err2 = overall_err_calc(RAc_err2, Dec_err2, corr2)


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
                "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                "ga15-odd-even-com-RAcerr.eps")
comparison_plot(Dec_err_com1, Dec_err_com2,
                "$\\sigma_\\delta$ (mas) - odd",
                "$\\sigma_\\delta$ (mas) - even",
                "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                "ga15-odd-even-com-Decerr.eps")
comparison_plot(num_ses_com1, num_ses_com2,
                "$N_{session}$ - odd",
                "$N_{session}$ - even",
                "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                "ga15-odd-even-com-Numses.eps")
comparison_plot(num_obs_com1, num_obs_com2,
                "$N_{observation}$ - odd",
                "$N_{observation}$ - even",
                "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                "ga15-odd-even-com-Numobs.eps")
comparison_plot(sig_pos_max_com1, sig_pos_max_com2,
                "Elliptical error (mas) - odd",
                "Elliptical error (mas) - even",
                "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                "ga15-odd-even-com-maxerr.eps")

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
fout = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
            "odd-even.cat_diff", "w")

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
     num_ses_meani, num_obs_meani) in zip(comsou, RA_com1, Dec_com1,
                                          dRAc, dDec, dRAc_err, dDec_err, num_ses_mean, num_obs_mean):
    print("%-8s  %+10.4f  %+10.4f  %+8.4f  %+8.4f  %8.4f  %8.4f  %6d  %8d" %
          (comsoui, RA_com1i, Dec_com1i, dRAci, dDeci, dRAc_erri, dDec_erri,
           num_ses_meani, num_obs_meani), file=fout)
fout.close()


# Plot for differences vs error
difference_vs_error_plot(dRAc_err_2, dRAc**2,
                         "$\\sigma_{\\Delta\\alpha^*}^2 (mas^2)$",
                         "$\\Delta^2_{\\alpha^*} (mas^2)$",
                         "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                         "ga15-odd-even-dRAc-dRAcerr.eps")
difference_vs_error_plot(dDec_err_2, dDec**2,
                         "$\\sigma_{\\Delta\\delta}^2 (mas^2)$",
                         "$\\Delta^2_{\\delta} (mas^2)$",
                         "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                         "ga15-odd-even-dDec-dDecerr.eps")
difference_vs_error_plot(dRAc_err_2, dRAc,
                         "$\\sigma_{\\Delta\\alpha^*}^2 (mas^2)$",
                         "$\\Delta^2_{\\alpha^*} (mas^2)$",
                         "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                         "ga15-odd-even-dRAc-dRAcerr.eps")
difference_vs_error_plot(dDec_err_2, dDec**2,
                         "$\\sigma_{\\Delta\\delta}^2 (mas^2)$",
                         "$\\Delta^2_{\\delta} (mas^2)$",
                         "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                         "ga15-odd-even-dDec-dDecerr.eps")

# --------------------------------- END --------------------------------
