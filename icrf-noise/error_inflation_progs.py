#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: error_inflation_progs.py
"""
Created on Tue Jun 26 12:38:35 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Functions to inflate the formal error, estimate the noise floor and scale factor.

"""

import numpy as np
from numpy import cos, deg2rad, sqrt
from scipy.optimize import curve_fit


__all__ = ["error_inflation", "wrms_calc", "sf_nf_calc", "nf_sf_calc_LSQ",
           "pos_offset_wrms_soubinned"]



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
    """Calculate the wrms of an array with removing the bias.
    """
    
    p = 1. / x_err**2
    p_tot = np.sum(p)
    mean = np.dot(x, p) / p_tot
    wrms = sqrt(np.sum((x-mean)**2 * p) / p_tot)
    
    return wrms


def sf_nf_calc(x, x_err):
    """Calculate the scale factor and noise floor.
    
    Parameters
    ----------
    x : array_like of float type
        estimate
    dx_err : array_like of float type
        formal error

    Returns
    ----------
    sf : float
        scaling fator
    nf : float
        noise floor    
    """
    
    # Calculate the scale factor
    sf = np.std(x / x_err)
#     m = np.mean(nor_dRAc[con])
    
    # Calculate the noise floor
    nf = wrms_calc(x, x_err) / sqrt(2)
    
    return sf, nf


def nf_sf_calc_LSQ(num_mean, dRAc, dDec, dRAc_err, dDec_err,
                      interv_size=50, min_num=100,
                      init_vals_RAc=[1.6, 0.2], init_vals_Dec=[1.6, 0.2]):
    """Estimate scaling factor and noise floor based on declination.
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

    # For Dec. component
    # init_vals_Dec = [1.6, 0.2]
    best_vals_Dec, covar_Dec = curve_fit(
        error_inflation, dDec_err_med_s, dDec_wrms_s, p0=init_vals_Dec)
    sf_Dec, nf_Dec = best_vals_Dec
    sf_Dec_err, nf_Dec_err = np.sqrt(np.diag(covar_Dec))
    # post-fit residual(O -C)
    res_Dec = dDec_wrms_s - error_inflation(
        dDec_err_med_s, sf_Dec, nf_Dec)
    chi2_Dec = np.dot(res_Dec, res_Dec) * 1.e6
    red_chi2_Dec = chi2_Dec / (dDec_err_med_s.size - 2)

    return [sf_RA, sf_RA_err, nf_RA, nf_RA_err,
            sf_Dec, sf_Dec_err, nf_Dec, nf_Dec_err]


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

        x_wrms[i] = wrms_calc(xs, x_errs) / sqrt(2.)
        err_med[i] = np.median(x_errs) / sqrt(2.)
        num_min[i] = np.median(num_sorts)

    return x_wrms, err_med, num_min

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
#         x_wrms[i] = np.std(xs)
        err_mid[i] = (interv_b + interv_e) / 2.
        err_med[i] = np.median(x_errs)
        
        print("(%.3f, %.3f] mas:  %4d sources"%(interv_b, interv_e, xs.size))

    return x_wrms, err_med, err_mid


def declination_error_calc_DSM(Dec, dRAc, dDec, dRAc_err, dDec_err):
    """Estimate decliantion-dependent scaling factor and noise floor.

    Parameters
    ----------
    Dec : array_like of float type
        declination in degree
    dRAc/dDec : array_like of float type
        offset in right ascension/declination component in mas
    dRAc_err/dDec_err : array_like of float type
        formal uncertainty of dRAc/dDec in mas

    Returns
    -------
    """

    # -90 ~ +90
    bin_size = 15
    dec_nodes = np.arange(-90, 90, bin_size)

    dec_means = dec_nodes + bin_size / 2.
    sf_RAs = np.zeros_like(dec_means)
    nf_RAs = np.zeros_like(dec_means)
    sf_Decs = np.zeros_like(dec_means)
    nf_Decs = np.zeros_like(dec_means)
    bin_nums = np.zeros_like(dec_means)

    print("#  bin   Num_sou    sf_ra    nf_ra   sf_dec   sf_dec")

    for i, dec_node in enumerate(dec_nodes):
        # Divide into sub-group
        con = (Dec >= dec_node) & (Dec < dec_node + bin_size)
        dRAc_s = dRAc[con]
        dDec_s = dDec[con]
        dRAc_err_s = dRAc_err[con]
        dDec_err_s = dDec_err[con]

        # Caculate scale factors and noise floor
        sf_RA, nf_RA = sf_nf_calc(dRAc_s, dRAc_err_s)
        sf_Dec, nf_Dec = sf_nf_calc(dDec_s, dDec_err_s)
        
        print("[%+3d,%+3d)  %4d  %7.3f  %7.3f  %7.3f  %7.3f" %
              (dec_node, dec_node + bin_size, dRAc_s.size,
               sf_RA, nf_RA, sf_Dec, nf_Dec))
        
        sf_RAs[i], nf_RAs[i] = sf_RA, nf_RA
        sf_Decs[i], nf_Decs[i] = sf_Dec, nf_Dec
        bin_nums[i] =  dRAc_s.size
        
    return dec_means, bin_nums, sf_RAs, nf_RAs, sf_Decs, nf_Decs
