#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: VLBI_internal_error.py
"""
Created on Thu May 17 11:31:14 2018

@author: Neo(liuniu@smail.nju.edu.cn)


This code is used to determine the scaling factor and noise floor of the
VLBI formal error.
"""

import numpy as np
from numpy import cos, deg2rad, sqrt
from scipy.optimize import curve_fit2
import time
# My modules
from list_crossmatch import list_crossmatch
from error_plot import (comparison_plot, comparison_plot2,
                        comparison_plot3, comparison_plot3_1,
                        difference_vs_error_plot, sf_nf_dec_plot)


# -----------------------------  FUNCTIONS -----------------------------
def error_inflation(dx_err, s, f):
    """dx_wrms^2 ~ (s * dx_err)^2 + f^2

    Parameters
    ----------
    dx_err : array_like of float type
        combined formal error
    s : float
        scaling fator
    f : float
        noise floor

    Returns
    ----------
    dx_wrms : scatter of positonal offset
    """

    return sqrt((s * dx_err)**2 + f**2)


def wrms_calc(x, x_err):
    """Calculate the wrms of a set of data.

    Parameters
    ----------
    x : array_like of float
        data
    x_err : array_like of float
        formal error of x

    Returns
    ----------
    wrms : float
        wrms
    """

    # return sqrt(np.dot(x/x_err, x/x_err) / np.dot(1./x_err, 1./x_err))
    return sqrt(np.sum(x**2/x_err**2) / np.sum(1./x_err**2))


def offset_wrms_errbinned(x, x_err, interv=0.04, max_err=1):
    """Calculate the scatter of positional offset binned by formal error

    Parameters
    ----------
    x : array_like of float
        data
    x_err : array_like of float
        formal error of x
    interv : float
        width of interval in mas
    max_err : float
        maximum of x_err in mas

    Returns
    ----------
    x_wrms : float of array
        wrms of x binned by interval of x_err
    err_med : float of array
        dedian of interval of x_err
    err_mid : float of array
        mid-point of interval of x_err
    """

    # interv = 40  # uas
    # max_err = 1000  # uas
    num_interv = int(max_err / interv)
    x_wrms = np.zeros(num_interv)
    err_mid = np.ones(num_interv)
    err_med = np.ones(num_interv)

    for i in range(num_interv):
        interv_b = 0.001 + i * interv
        interv_e = interv_b + interv
        con = (x_err > interv_b) & (x_err <= interv_e)
        xs = x[con]
        x_errs = x_err[con]

        x_wrms[i] = wrms_calc(xs, x_errs)
        err_mid[i] = (interv_b + interv_e) / 2.
        err_med[i] = np.median(x_errs)

    return x_wrms, err_med, err_mid


def pos_offset_wrms_soubinned(num_sort, x, x_err, interv_size=50):
    """Calculate the scatter of positional offset

    Parameters
    ----------
    num_sort : array_like of int
        number of sessions/observations a source was observed
    x : array_like of float
        data
    x_err : array_like of float
        formal error of x
    interv_size : int
        number of source in an interval

    Returns
    ----------
    x_wrms : array of float
        wrms of x binned by interval of x_err
    err_med : array of float
        mdedian of x_err in every interval
    num_min : array of int
        minimum number of session/observation in the subset
    """

    if x.size <= interv_size:
        print("# Too smaller sample!")
        exit()
    else:
        interv_num = x.size - interv_size + 1

    x_wrms = np.zeros(interv_num)
    err_med = np.ones(interv_num)
    num_min = np.ones(interv_num)

    for i in range(interv_num):
        ind_b, ind_e = i, i + interv_size
        num_sorts = num_sort[ind_b: ind_e]
        xs = x[ind_b: ind_e]
        x_errs = x_err[ind_b: ind_e]

        x_wrms[i] = wrms_calc(xs, x_errs)
        err_med[i] = np.median(x_errs)
        num_min[i] = min(num_sorts)

    return x_wrms, err_med, num_min


def sorted_error_calc(num_mean, dRAc, dDec, dRAc_err, dDec_err, flog,
                      interv_size=50, min_num=100,
                      init_vals_RAc=[1.6, 0.2], init_vals_Dec=[1.6, 0.2],
                      fig_out=""):
    """Estimate scaling factor and noise floor based on session or obs. number.
    """

    ind = np.argsort(num_mean)
    num_sort = np.take(num_mean, ind)
    dRAc_sort = np.take(dRAc, ind)
    dDec_sort = np.take(dDec, ind)
    dRAc_err_sort = np.take(dRAc_err, ind)
    dDec_err_sort = np.take(dDec_err, ind)

# Calculate the wrms of positional offset
    dRAc_wrms, dRAc_err_med, ses_min_RA = pos_offset_wrms_soubinned(
        num_sort, dRAc_sort, dRAc_err_sort, interv_size)
    dDec_wrms, dDec_err_med, ses_min_Dec = pos_offset_wrms_soubinned(
        num_sort, dDec_sort, dDec_err_sort, interv_size)

    # Applied condition on number of session
    # # use a limit on N_ses
    # min_num = 100
    con = (ses_min_RA > min_num)
    # con = (ses_min_RA >= 10) & (ses_min_RA <= 20)
    dRAc_err_med_s = dRAc_err_med[con]
    dRAc_wrms_s = dRAc_wrms[con]
    dDec_err_med_s = dDec_err_med[con]
    dDec_wrms_s = dDec_wrms[con]

    # init_vals = [1.6, 0.2]
    best_vals_RA, covar_RA = curve_fit(
        error_inflation, dRAc_err_med_s, dRAc_wrms_s, p0=init_vals_RAc)
    sf_RA, nf_RA = best_vals_RA
    sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))
    # post-fit residual(O -C)
    res_RA = dRAc_wrms_s - error_inflation(dRAc_err_med_s, sf_RA, nf_RA)
    chi2_RA = np.dot(res_RA, res_RA) * 1.e6
    red_chi2_RA = chi2_RA / (dRAc_err_med_s.size - 2)

    print("# For R.A. component (N >= %d):\n"
          "# scaling factor: %.2f +/- %.2f\n"
          "# noise floor: %.3f +/- %.3f mas\n"
          "# data points: %d\n"
          "# variables: 2\n"
          "# chi-square: %.2f (uas^2)\n"
          "# reduced chi-square: %.2f (uas^2)\n" %
          (min_num, sf_RA, sf_RA_err, nf_RA, nf_RA_err,
           dRAc_err_med_s.size, chi2_RA, red_chi2_RA), file=flog)

    # For Dec. component
    # init_vals_Dec = [1.6, 0.2]
    best_vals_Dec, covar_Dec = curve_fit(
        error_inflation, dDec_err_med_s, dDec_wrms_s, p0=init_vals_Dec)
    sf_Dec, nf_Dec = best_vals_Dec
    sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))
    # post-fit residual(O -C)
    res_Dec = dDec_wrms_s - error_inflation(
        dDec_err_med_s, sf_Dec, nf_Dec/1.e3)
    chi2_Dec = np.dot(res_Dec, res_Dec) * 1.e6
    red_chi2_Dec = chi2_Dec / (dDec_err_med_s.size - 2)
    print("For Dec. component (N >= %d):\n"
          "scaling factor: %.2f +/- %.2f\n"
          "noise floor: %.3f +/- %.3f mas\n"
          "# data points: %d\n"
          "# variables: 2\n"
          "# chi-square: %.2f (uas^2)\n"
          "# reduced chi-square: %.2f (uas^2)\n" %
          (min_num, sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err,
           dDec_err_med_s.size, chi2_Dec, red_chi2_Dec), file=flog)

    comparison_plot3_1(ses_min_RA, dRAc_wrms, dRAc_err_med,
                       "$N_{session}$",
                       "WRMS in $\\alpha\\,\\cos\\delta$ (mas)",
                       "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                       "%s-RA.eps" % fig_out, s=sf_RA, f=nf_RA)
    # "%s-RA.eps" % fig_out, s=sf_RA, f=nf_RA, ylog=True)
    comparison_plot3_1(ses_min_Dec, dDec_wrms, dDec_err_med,
                       "$N_{session}$",
                       "WRMS in $\\delta$(mas)",
                       "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                       "%s-Dec.eps" % fig_out, s=sf_Dec, f=nf_Dec)
    # "%s-Dec.eps" % fig_out, s=sf_Dec, f=nf_Dec, ylog=True)

    return [sf_RA, sf_RA_err, nf_RA, nf_RA_err,
            sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err]


# def obs_sorted_calc(num_ses_sort, dRAc_sort, dRAc_err_sort,
#                     interv_size=50, min_ses=100):
#     """Estimate scaling factor and noise floor based on session number.
#     """


# --------------------------------- MAIN -------------------------------
# Load data
datfile = ("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
           "odd-even.cat_diff")
sou = np.genfromtxt(datfile, usecols=(0,), dtype=str)
RA, Dec, dRAc, dDec, dRAc_err, dDec_err = np.genfromtxt(
    datfile, usecols=range(1, 7), unpack=True)
num_ses_mean, num_obs_mean = np.genfromtxt(
    datfile, dtype=int, usecols=range(7, 9), unpack=True)

# estimation of scaling factor and noise floor
# Method #01)
print("## 01) direct estimation:")

# For R.A. component
init_vals = [1.0, 0.00]
best_vals_RA, covar_RA = curve_fit(
    error_inflation, dRAc_err_2, dRAc**2, p0=init_vals)
sf_RA, nf_RA = best_vals_RA
sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))
print("For R.A. component:\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_RA, sf_RA_err, nf_RA, nf_RA_err))

# use a limit of err < 1 mas
con = (dRAc_err_2 <= 1) & (dRAc_err_2 >= 1.e-6)
dRAc_err_s = dRAc_err[con]
dRAc_s = np.fabs(dRAc[con])
best_vals_RA, covar_RA = curve_fit(
    error_inflation, dRAc_err_s, dRAc_s, p0=init_vals)
sf_RA, nf_RA = best_vals_RA
sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))
print("For R.A. component (0.001 mas <= err <= 1 mas):\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_RA, sf_RA_err, nf_RA, nf_RA_err))

# For Dec. component
best_vals_Dec, covar_Dec = curve_fit(
    error_inflation, dDec_err_2, dDec, p0=init_vals)
sf_Dec, nf_Dec = best_vals_Dec
sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))
print("For Dec. component:\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err))

# use a limit of 0.001 mas <= err <= 1 mas
con = (dDec_err_2 <= 1) & (dDec_err_2 >= 1.e-6)
dDec_err_s = dDec_err[con]
dDec_s = np.fabs(dDec[con])
best_vals_Dec, covar_Dec = curve_fit(
    error_inflation, dDec_err_s, dDec_s, p0=init_vals)
sf_Dec, nf_Dec = best_vals_Dec
sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))
print("For Dec. component (0.001 mas <= dDec_err <= 1 mas):\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err))


# # Method #02) Lambert AA 2014
print("## 02) Lambert AA 2014:")
# Calculate the positional offset scatter
# interv = 0.02
dRAc_wrms, dRAc_err_med, dRAc_err_mid = offset_wrms_errbinned(
    dRAc, dRAc_err, interv=0.01, max_err=1)
dDec_wrms, dDec_err_med, dDec_err_mid = offset_wrms_errbinned(
    dDec, dDec_err, interv=0.01, max_err=1)


# For R.A. component
best_vals_RA, covar_RA = curve_fit(
    error_inflation, dRAc_err_med, dRAc_wrms, p0=init_vals)
sf_RA, nf_RA = best_vals_RA
sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))

con = (dRAc_err_med <= 0.1)
dRAc_err_med_s = dRAc_err_med[con]
dRAc_wrms_s = dRAc_wrms[con]

best_vals_RA, covar_RA = curve_fit(
    error_inflation, dRAc_err_med_s, dRAc_wrms_s, p0=init_vals)
sf_RA, nf_RA = best_vals_RA
sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))

# # post-fit residual
# rsd_RA = dRAc_wrms - error_inflation(dRAc_err_med, sf_RA, nf_RA)
# rsd_std_RA = np.std(rsd_RA)

# # 3-sigma priciple
# con = (np.fabs(rsd_RA) <= 3 * rsd_std_RA) & (dRAc_err_med <= 0.1)
# dRAc_err_med_s = dRAc_err_med[con]
# dRAc_wrms_s = dRAc_wrms[con]

# best_vals_RA, covar_RA = curve_fit(
#     error_inflation, dRAc_err_med_s, dRAc_wrms_s, p0=init_vals)

# sf_RA, nf_RA = best_vals_RA
# sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))
dRAc_err_inf = error_inflation(dRAc_err_med, sf_RA, nf_RA)

print("For R.A. component:\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_RA, sf_RA_err, nf_RA, nf_RA_err))

# For Dec. component
# print(dDec_wrms - dDec_err_med)
# 8

best_vals_Dec, covar_Dec = curve_fit(
    error_inflation, dDec_err_med, dDec_wrms, p0=init_vals)
sf_Dec, nf_Dec = best_vals_Dec
sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))

con = (dDec_err_med <= 0.07)
dDec_err_med_s = dDec_err_med[con]
dDec_wrms_s = dDec_wrms[con]
# Manually remove 8th data point
# dDec_err_med_s[4] = dDec_wrms_s[4]
# dDec_err_med_s[7] = dDec_wrms_s[7]

best_vals_Dec, covar_Dec = curve_fit(
    error_inflation, dDec_err_med_s, dDec_wrms_s, p0=init_vals)
sf_Dec, nf_Dec = best_vals_Dec
sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))

# # post-fit residual
# rsd_Dec = dDec_wrms - error_inflation(dDec_err_med, sf_Dec, nf_Dec)
# rsd_std_Dec = np.std(rsd_Dec)

# # 3-sigma priciple
# con = (np.fabs(rsd_Dec) <= 3 * rsd_std_Dec) & (dDec_err_med <= 0.1)
# dDec_err_med_s = dDec_err_med[con]
# dDec_wrms_s = dDec_wrms[con]

# best_vals_Dec, covar_Dec = curve_fit(
#     error_inflation, dDec_err_med_s, dDec_wrms_s, p0=init_vals)
# sf_Dec, nf_Dec = best_vals_Dec
# sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))
dDec_err_inf = error_inflation(dDec_err_med, sf_Dec, nf_Dec)

print("For Dec. component:\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err))


# Plot
comparison_plot2(dRAc_err_med, dRAc_wrms, dRAc_err_inf,
                 dDec_err_med, dDec_wrms, dDec_err_inf,
                 "Formal error in $\\alpha\\,\\cos\\delta$ (mas)",
                 # "Scatter of the offset in $\\alpha\\,\\cos\\delta$ (mas)",
                 "Scatter in $\\alpha\\,\\cos\\delta$ (mas)",
                 "Formal error in $\\delta$ (mas)",
                 # "Scatter of the offset in $\\delta$ (mas)",
                 "Scatter in $\\delta$ (mas)",
                 # "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                 # "ga15-odd-even-com-wrms-combinederr.eps")
                 "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                 "ga15-odd-even-com-wrms-combinederr-log.eps")

"""
## 02) Lambert AA 2014:
For R.A. component:
scaling factor: 1.75 +/- 0.14
noise floor: 0.005 +/- 0.060 mas

For Dec. component:
scaling factor: 1.53 +/- 0.19
noise floor: 0.019 +/- 0.016 mas
"""

# # Method #03) IERS TN 35/ ICRF2
print("## 03) IERS TN 35/ ICRF2:")

# # Sort the data according to the number of sessions/observations
print("##     session-sorted:")
ind = np.argsort(num_ses_mean)
num_ses_sort = np.take(num_ses_mean, ind)
num_obs_sort = np.take(num_obs_mean, ind)
dRAc_sort = np.take(dRAc, ind)
dDec_sort = np.take(dDec, ind)
dRAc_err_sort = np.take(dRAc_err, ind)
dDec_err_sort = np.take(dDec_err, ind)

# Calculate the positional offset scatter
dRAc_wrms, dRAc_err_med, ses_min_RA = pos_offset_wrms_soubinned(
    num_ses_sort, dRAc_sort, dRAc_err_sort, interv_size=50)
dDec_wrms, dDec_err_med, ses_min_Dec = pos_offset_wrms_soubinned(
    num_ses_sort, dDec_sort, dDec_err_sort, interv_size=50)

# session-sorted
comparison_plot3_1(ses_min_RA, dRAc_wrms, dRAc_err_med,
                   "$N_{session}$",
                   "WRMS in $\\alpha\\,\\cos\\delta$ (mas)",
                   "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                   "ga15-odd-even-com-wrms-combinederr-sortses-RA.eps")
comparison_plot3_1(ses_min_Dec, dDec_wrms, dDec_err_med,
                   "$N_{session}$",
                   "WRMS in $\\delta$(mas)",
                   "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                   "ga15-odd-even-com-wrms-combinederr-sortses-Dec.eps")

# For R.A. component
best_vals_RA, covar_RA = curve_fit(
    error_inflation, dRAc_err_med, dRAc_wrms, p0=init_vals)
sf_RA, nf_RA = best_vals_RA
sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))
print("For R.A. component:\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_RA, sf_RA_err, nf_RA, nf_RA_err))

# For Dec. component
best_vals_Dec, covar_Dec = curve_fit(
    error_inflation, dDec_err_med, dDec_wrms, p0=init_vals)
sf_Dec, nf_Dec = best_vals_Dec
sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))
print("For Dec. component:\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n" %
      (sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err))


# Calculate the wrms of positional offset
dRAc_wrms, dRAc_err_med, ses_min_RA = pos_offset_wrms_soubinned(
    num_ses_sort, dRAc_sort, dRAc_err_sort, interv_size=50)
dDec_wrms, dDec_err_med, ses_min_Dec = pos_offset_wrms_soubinned(
    num_ses_sort, dDec_sort, dDec_err_sort, interv_size=50)

# Applied condition on number of session
# # use a limit on N_ses
min_ses = 100
con = (ses_min_RA > min_ses)
# con = (ses_min_RA >= 10) & (ses_min_RA <= 20)
dRAc_err_med_s = dRAc_err_med[con]
dRAc_wrms_s = dRAc_wrms[con]
dDec_err_med_s = dDec_err_med[con]
dDec_wrms_s = dDec_wrms[con]

init_vals = [1.6, 0.2]
best_vals_RA, covar_RA = curve_fit(
    error_inflation, dRAc_err_med_s, dRAc_wrms_s, p0=init_vals)
sf_RA, nf_RA = best_vals_RA
sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))
# post-fit residual(O -C)
res_RA = dRAc_wrms_s - error_inflation(dRAc_err_med_s, sf_RA, nf_RA)
chi2_RA = np.dot(res_RA, res_RA) * 1.e6
red_chi2_RA = chi2_RA / (dRAc_err_med_s.size - 2)

print("# For R.A. component (N_ses >= %d):\n"
      "# scaling factor: %.2f +/- %.2f\n"
      "# noise floor: %.3f +/- %.3f mas\n"
      "# data points: %d\n"
      "# variables: 2\n"
      "# chi-square: %.2f (uas^2)\n"
      "# reduced chi-square: %.2f (uas^2)\n" %
      (min_ses, sf_RA, sf_RA_err, nf_RA, nf_RA_err,
       dRAc_err_med_s.size, chi2_RA, red_chi2_RA))

# For Dec. component
init_vals = [1.6, 0.2]
best_vals_Dec, covar_Dec = curve_fit(
    error_inflation, dDec_err_med_s, dDec_wrms_s, p0=init_vals)
sf_Dec, nf_Dec = best_vals_Dec
sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))

# post-fit residual(O -C)
res_Dec = dDec_wrms_s - error_inflation(dDec_err_med_s, sf_Dec, nf_Dec/1.e3)
chi2_Dec = np.dot(res_Dec, res_Dec) * 1.e6
red_chi2_Dec = chi2_Dec / (dDec_err_med_s.size - 2)
print("For Dec. component (N_ses >= %d):\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n"
      "# data points: %d\n"
      "# variables: 2\n"
      "# chi-square: %.2f (uas^2)\n"
      "# reduced chi-square: %.2f (uas^2)\n" %
      (min_ses, sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err,
       dDec_err_med_s.size, chi2_Dec, red_chi2_Dec))

# Plot
comparison_plot3_1(ses_min_RA, dRAc_wrms, dRAc_err_med,
                   "$N_{session}$",
                   "WRMS in $\\alpha\\,\\cos\\delta$ (mas)",
                   "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                   # "ga15-odd-even-com-wrms-combinederr-sortses-RA.eps",
                   "ga15-odd-even-com-wrms-combinederr-sortses-RA-log.eps",
                   s=sf_RA, f=nf_RA)
comparison_plot3_1(ses_min_Dec, dDec_wrms, dDec_err_med,
                   "$N_{session}$",
                   "WRMS in $\\delta$(mas)",
                   "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                   # "ga15-odd-even-com-wrms-combinederr-sortses-Dec.eps",
                   "ga15-odd-even-com-wrms-combinederr-sortses-Dec-log.eps",
                   s=sf_Dec, f=nf_Dec)


print("##     observstion-sorted:")
ind = np.argsort(num_obs_mean)
num_ses_sort = np.take(num_ses_mean, ind)
num_obs_sort = np.take(num_obs_mean, ind)
dRAc_sort = np.take(dRAc, ind)
dDec_sort = np.take(dDec, ind)
dRAc_err_sort = np.take(dRAc_err, ind)
dDec_err_sort = np.take(dDec_err, ind)

# Calculate the wrms of positional offset
dRAc_wrms, dRAc_err_med, obs_min_RA = pos_offset_wrms_soubinned(
    num_obs_sort, dRAc_sort, dRAc_err_sort, interv_size=50)
dDec_wrms, dDec_err_med, obs_min_Dec = pos_offset_wrms_soubinned(
    num_obs_sort, dDec_sort, dDec_err_sort, interv_size=50)

# use a limit on N_obs
min_obs = 1000
con = (obs_min_RA > min_obs)
# con = (obs_min_RA >= 10) & (obs_min_RA <= 20)
dRAc_err_med_s = dRAc_err_med[con]
dRAc_wrms_s = dRAc_wrms[con]
dDec_err_med_s = dDec_err_med[con]
dDec_wrms_s = dDec_wrms[con]

# For RA component
init_vals = [1.6, 0.0]
best_vals_RA, covar_RA = curve_fit(
    error_inflation, dRAc_err_med_s, dRAc_wrms_s, p0=init_vals)
sf_RA, nf_RA = best_vals_RA
sf_RA_err, nf_RA_err = np.sqrt(np.diag(covar_RA))
# post-fit residual(O -C)
res_RA = dRAc_wrms_s - error_inflation(dRAc_err_med_s, sf_RA, nf_RA)
chi2_RA = np.dot(res_RA, res_RA) * 1.e6
red_chi2_RA = chi2_RA / (dRAc_err_med_s.size - 2)

print("# For R.A. component (N_obs >= %d):\n"
      "# scaling factor: %.2f +/- %.2f\n"
      "# noise floor: %.3f +/- %.3f mas\n"
      "# data points: %d\n"
      "# variables: 2\n"
      "# chi-square: %.2f (uas^2)\n"
      "# reduced chi-square: %.2f (uas^2)\n" %
      (min_obs, sf_RA, sf_RA_err, nf_RA, nf_RA_err,
       dRAc_err_med_s.size, chi2_RA, red_chi2_RA))

# For Dec. component
init_vals = [1.6, 1.]
best_vals_Dec, covar_Dec = curve_fit(
    error_inflation, dDec_err_med_s, dDec_wrms_s, p0=init_vals)
sf_Dec, nf_Dec = best_vals_Dec
sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))
# post-fit residual(O -C)
res_Dec = dDec_wrms_s - error_inflation(dDec_err_med_s, sf_Dec, nf_Dec/1.e3)
chi2_Dec = np.dot(res_Dec, res_Dec) * 1.e6
red_chi2_Dec = chi2_Dec / (dDec_err_med_s.size - 2)
print("For Dec. component (N_obs >= %d):\n"
      "scaling factor: %.2f +/- %.2f\n"
      "noise floor: %.3f +/- %.3f mas\n"
      "# data points: %d\n"
      "# variables: 2\n"
      "# chi-square: %.2f (uas^2)\n"
      "# reduced chi-square: %.2f (uas^2)\n" %
      (min_obs, sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err,
       dDec_err_med_s.size, chi2_Dec, red_chi2_Dec))

# Plot
comparison_plot3_1(obs_min_RA, dRAc_wrms, dRAc_err_med,
                   "$N_{observation}$",
                   "WRMS in $\\alpha\\,\\cos\\delta$ (mas)",
                   "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                   "ga15-odd-even-com-wrms-combinederr-sortobs-RA.eps",
                   s=sf_RA, f=nf_RA)
# "ga15-odd-even-com-wrms-combinederr-sortobs-RA-log.eps",
# s=sf_Dec, f=nf_Dec, ylog=True)

comparison_plot3_1(obs_min_Dec, dDec_wrms, dDec_err_med,
                   "$N_{observation}$",
                   "WRMS in $\\delta$(mas)",
                   "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                   "ga15-odd-even-com-wrms-combinederr-sortobs-Dec.eps",
                   s=sf_Dec, f=nf_Dec)
# "ga15-odd-even-com-wrms-combinederr-sortobs-Dec-log.eps",
# s=sf_Dec, f=nf_Dec, ylog=True)


# Assess the declination-dependent scaling factor and noise floor
flog = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/logs/"
            "VLBI_formalerror_dec.log", "w")
fout = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/logs/"
            "VLBI_formalerror_dec.dat", "w")


# -90 ~ +90
bin_size = 15
interv_sizes = np.array([2, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10])
min_nums = np.array([5, 20, 20, 40, 100, 100, 50, 50, 50, 20, 30, 2])
dec_nodes = np.arange(-90, 90, bin_size)

print("# Declination(Deg)              R.A.              "
      "       Dec.\n"
      "#                  -----------------------------  "
      "-----------------------------\n"
      "#                  ScaleFactor NoiseFloor(mas)    "
      "ScaleFactor NoiseFloor(mas)")

for (dec_node, interv_size, min_num) in zip(dec_nodes, interv_sizes, min_nums):
    con = (Dec_com2 >= dec_node) & (Dec_com2 < dec_node + bin_size)
    print("# Number of source in [%+3d, %+3d): %4d\n" %
          (dec_node, dec_node + bin_size, Dec_com2[con].size),
          file=flog)
    dec_mean = dec_node + bin_size / 2.

    # Divide into sub-group
    dRAc_s = dRAc[con]
    dDec_s = dDec[con]
    dRAc_err_s = dRAc_err[con]
    dDec_err_s = dDec_err[con]
    num_ses_mean_s = num_ses_mean[con]
    num_obs_mean_s = num_obs_mean[con]

    [sf_RA, sf_RA_err, nf_RA, nf_RA_err,
     sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err] = sorted_error_calc(
        num_ses_mean_s, dRAc_s, dDec_s,
        dRAc_err_s, dDec_err_s, flog,
        interv_size, min_num,
        init_vals_RAc=[1.0, 0.0], init_vals_Dec=[1.0, 0.0],
        fig_out="%+d_%+d" % (dec_node, dec_node + bin_size))

    print("%4d  %+5.1f  %5.3f  %6.4f  %5.3f  %6.4f" %
          (dRAc_s.size, dec_mean, sf_RA, nf_RA, sf_Dec, nf_Dec), file=fout)

    print("[%+3d, %+3d)"
          "  %+5.3f +/- %5.3f  %+5.4f +/- %7.4f"
          "  %+5.3f +/- %5.3f  %+5.4f +/- %7.4f" %
          (dec_node, dec_node + bin_size,
           sf_RA, sf_RA_err, nf_RA, nf_RA_err,
           sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err))
flog.close()
fout.close()
"""
# Declination(Deg)              R.A.                   Dec.
#                  -------------------------  ---------------------------------
#                  ScaleFactor NoiseFloor(mas)  ScaleFactor NoiseFloor(mas)
[-90, -75)  +1.368 +/- 0.483  +0.0938 +/-  0.1000  +1.454 +/- 0.325  +0.1007 +/-  0.0742
[-75, -60)  +1.050 +/- 0.166  +0.0175 +/-  0.0339  +1.206 +/- 0.083  +0.0000 +/- 60.8590
[-60, -45)  +1.088 +/- 0.053  +0.0215 +/-  0.0083  +1.451 +/- 0.065  +0.0221 +/-  0.0300
[-45, -30)  +1.174 +/- 0.045  +0.0217 +/-  0.0025  +0.586 +/- 0.122  +0.0421 +/-  0.0046
[-30, -15)  +1.399 +/- 0.074  -0.0000 +/- 13.9153  +1.139 +/- 0.063  +0.0180 +/-  0.0025
[-15,  +0)  +1.145 +/- 0.068  +0.0000 +/- 10.4830  +0.797 +/- 0.123  +0.0210 +/-  0.0018
[ +0, +15)  +1.583 +/- 0.041  +0.0034 +/-  0.0072  +0.990 +/- 0.090  +0.0304 +/-  0.0030
[+15, +30)  +1.438 +/- 0.050  +0.0134 +/-  0.0020  +0.980 +/- 0.046  +0.0130 +/-  0.0022
[+30, +45)  +1.088 +/- 0.098  +0.0084 +/-  0.0026  +0.610 +/- 0.230  +0.0265 +/-  0.0021
[+45, +60)  +1.899 +/- 0.118  +0.0000 +/- 14.9117  +1.406 +/- 0.040  +0.0201 +/-  0.0016
[+60, +75)  +1.101 +/- 0.057  +0.0163 +/-  0.0016  +2.025 +/- 0.063  +0.0004 +/-  0.1243
[+75, +90)  +1.121 +/- 0.206  +0.0305 +/-  0.0095  +0.759 +/- 0.087  +0.0193 +/-  0.0046
"""


dec, sf_RA, nf_RA, sf_Dec, nf_Dec = np.genfromtxt(
    "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/logs/"
    "VLBI_formalerror_dec.dat", usecols=range(1, 6), unpack=True)
sf_nf_dec_plot(dec, sf_RA, nf_RA, sf_Dec, nf_Dec,
               "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
               "ga15-odd-even-nf-sf-dec.eps")
# --------------------------------- END --------------------------------
