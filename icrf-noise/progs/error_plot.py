#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: error_plot.py
"""
Created on Sat May 19 19:48:23 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import matplotlib.pyplot as plt


__all__ = ["error_vs_numses", "error_vs_numses2",
           "error_vs_numobs", "error_vs_numobs2",
           "maxerror_vs_numses", "maxerror_vs_numses2",
           "maxerror_vs_numobs", "maxerror_vs_numobs2",
           "overallerror_vs_numses", "overallerror_vs_numobs",
           "overallerror_vs_numses2", "overallerror_vs_numobs2",
           "comparison_plot", "difference_vs_error_plot"]


# -----------------------------  FUNCTIONS -----------------------------
def error_vs_numses(RAc_err, Dec_err, num_ses, output):
    """Plot of formal error vs number of sessions

    Parameters
    ----------
    RAc_err/Dec_err : array of float
        formal uncertainty of RA*cos(Dec)/Dec in the VLBI solution
    num_ses : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.plot(num_ses, RAc_err, ".", markersize=1)
    ax1.plot(num_ses, Dec_err, ".", markersize=1)

    N = np.arange(1, 1e4)
    sigma = 1. / np.sqrt(N)
    ax0.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")
    ax1.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax0.set_ylabel("$\\sigma_{\\alpha^*}$ (mas)")
    ax0.set_yscale("log")
    ax0.grid(True)

    ax1.set_ylabel("$\\sigma_\\delta$ (mas)")
    # ax1.set_ylim([-5, 50])
    ax1.set_xlabel("$N_{session}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def error_vs_numses2(RAc_err1, Dec_err1, num_ses1,
                     RAc_err2, Dec_err2, num_ses2, output):
    """Plot of formal error vs number of sessions(two catalogs)

    Parameters
    ----------
    RAc_err/Dec_err : array of float
        formal uncertainty of RA*cos(Dec)/Dec in the VLBI solution
    num_ses : array of int
        number of sessions used in the solution
    dat_dir : string
        directory of the data
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.plot(num_ses1, RAc_err1, "r.", markersize=1, label="odd")
    ax0.plot(num_ses2, RAc_err2, "bx", markersize=1, label="even")
    ax1.plot(num_ses1, Dec_err1, "r.", markersize=1, label="odd")
    ax1.plot(num_ses2, Dec_err2, "bx", markersize=1, label="even")

    N = np.arange(1, 1e4)
    sigma = 1. / np.sqrt(N)
    ax0.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")
    ax1.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax0.set_ylabel("$\\sigma_ {\\alpha^*}$ (mas)")
    ax0.set_yscale("log")
    ax0.grid(True)
    ax0.legend()

    ax1.set_ylabel("$\\sigma_\\delta$ (mas)")
    # ax1.set_ylim([-5, 50])
    ax1.set_xlabel("$N_{session}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(True)
    ax1.legend()

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def error_vs_numobs(RAc_err, Dec_err, num_obs, output):
    """Plot of formal error vs number of observations

    Parameters
    ----------
    RAc_err/Dec_err : array of float
        formal uncertainty of RA*cos(Dec)/Dec in the VLBI solution
    num_obs : array of int
        number of observations used in the solution
    output : string
        output filename, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.plot(num_obs, RAc_err, "r.", markersize=1)
    ax1.plot(num_obs, Dec_err, "r.", markersize=1)

    N = np.arange(1, 1e6)
    sigma = 1. / np.sqrt(N)
    ax0.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")
    ax1.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax0.set_ylabel("$\\sigma_{\\alpha^*}$ (mas)")
    ax0.set_yscale("log")
    ax0.grid(True)

    ax1.set_ylabel("$\\sigma_\\delta$ (mas)")
    # ax1.set_ylim([-5, 50])
    ax1.set_xlabel("$N_{observation}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def error_vs_numobs2(RAc_err1, Dec_err1, num_obs1,
                     RAc_err2, Dec_err2, num_obs2, output):
    """Plot of formal error vs number of observations(two catalogs)

    Parameters
    ----------
    RAc_err/Dec_err : array of float
        formal uncertainty of RA*cos(Dec)/Dec in the VLBI solution
    num_obs : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    # RA component
    ax0.plot(num_obs1, RAc_err1, "r.", markersize=1, label="odd")
    ax0.plot(num_obs2, RAc_err2, "bx", markersize=1, label="even")

    # offset_RAc = np.mean(np.log10(RAc_err1) + 0.5 * np.log10(num_obs1))
    # # print(offset_RAc)
    # """0.27008259242888144
    # """
    # offset_RAc = 0.7
    N = np.arange(1, 1e6)
    # sigma_RA = 10**(np.log10(offset_RAc) - 0.5 * np.log10(N))
    # ax0.plot(N, sigma_RA, "k", label="$y=5.0/\\sqrt{N}$")
    sigma = 1. / np.sqrt(N)
    ax0.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax0.set_ylabel("$\\sigma_{\\alpha^*}$ (mas)")
    ax0.set_yscale("log")
    ax0.grid(True)
    ax0.legend()

    ax1.plot(num_obs1, Dec_err1, "r.", markersize=1, label="odd")
    ax1.plot(num_obs2, Dec_err2, "bx", markersize=1, label="even")

    # offset_Dec = np.mean(np.log10(Dec_err1) + 0.5 * np.log10(num_obs1))
    # print(offset_Dec)
    """0.46463189493215595
    """
    # offset_Dec = 0.7
    # sigma_Dec = 10**(np.log10(offset_Dec) - 0.5 * np.log10(N))
    # ax1.plot(N, sigma_Dec, "k", label="$y=5.0/\\sqrt{x}$")
    ax1.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax1.set_ylabel("$\\sigma_\\delta$ (mas)")
    # ax1.set_ylim([-5, 50])
    ax1.set_xlabel("$N_{observation}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(True)
    ax1.legend()

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def maxerror_vs_numses(sig_pos_max, num_ses, output):
    """Plot of ellipical error vs number of sessions

    Parameters
    ----------
    sig_pos_max : array of float
        semi-major axis of error ellipse in the VLBI solution
    num_ses : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_ses, sig_pos_max, ".", markersize=1)

    N = np.arange(1, 1e4)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Semi-major axis of dispersion ellipse (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{session}$")
    ax.set_xscale("log")

    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def maxerror_vs_numobs(sig_pos_max, num_obs, output):
    """Plot of ellipical error vs number of observations

    Parameters
    ----------
    sig_pos_max : array of float
        semi-major axis of error ellipse in the VLBI solution
    num_obs : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_obs, sig_pos_max, ".", markersize=1)

    N = np.arange(1, 1e6)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Semi-major axis of dispersion ellipse (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{observation}$")
    ax.set_xscale("log")

    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def overallerror_vs_numses(overall_err, num_ses, output):
    """Plot of overall formal error vs number of sessions

    Parameters
    ----------
    overall_err : array of float
        overall formal error sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)
        in the VLBI solution
    num_ses : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_ses, overall_err, ".", markersize=1)

    N = np.arange(1, 1e4)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Overall formal error (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{session}$")
    ax.set_xscale("log")

    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def overallerror_vs_numobs(overall_err, num_obs, output):
    """Plot of ellipical error vs number of sessions

    Parameters
    ----------
    overall_err : array of float
        overall formal error sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)
        in the VLBI solution
    num_obs : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_obs, overall_err, ".", markersize=1)

    N = np.arange(1, 1e6)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Overall formal error (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{observation}$")
    ax.set_xscale("log")

    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def maxerror_vs_numses2(sig_pos_max1, num_ses1, sig_pos_max2, num_ses2,
                        output):
    """Plot of ellipical error vs number of sessions(2 catalog)

    Parameters
    ----------
    sig_pos_max : array of float
        semi-major axis of error ellipse in the VLBI solution
    num_ses : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_ses1, sig_pos_max1, "r.", markersize=1, label="odd")
    ax.plot(num_ses2, sig_pos_max2, "bx", markersize=1, label="even")

    N = np.arange(1, 1e4)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Semi-major axis of dispersion ellipse (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{session}$")
    ax.set_xscale("log")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def maxerror_vs_numobs2(sig_pos_max1, num_obs1, sig_pos_max2, num_obs2,
                        output):
    """Plot of ellipical error vs number of observations(two catalogs)

    Parameters
    ----------
    sig_pos_max : array of float
        semi-major axis of error ellipse in the VLBI solution
    num_obs : array of int
        number of obseravtions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_obs1, sig_pos_max1, "r.", markersize=1, label="odd")
    ax.plot(num_obs2, sig_pos_max2, "bx", markersize=1, label="even")

    N = np.arange(1, 1e6)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Semi-major axis of dispersion ellipse (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{observation}$")
    ax.set_xscale("log")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def overallerror_vs_numses2(overall_err1, num_ses1, overall_err2, num_ses2,
                            output):
    """Plot of overall formal error  vs number of sessions(2 catalogs)

    Parameters
    ----------
    overall_err : array of float
        overall formal error ellipse in the VLBI solution
    num_ses : array of int
        number of sessions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_ses1, overall_err1, "r.", markersize=1, label="odd")
    ax.plot(num_ses2, overall_err2, "bx", markersize=1, label="even")

    N = np.arange(1, 1e4)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Overall formal error (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{session}$")
    ax.set_xscale("log")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def overallerror_vs_numobs2(overall_err1, num_obs1, overall_err2, num_obs2,
                            output):
    """Plot of overall formal error vs number of observations(two catalogs)

    Parameters
    ----------
    overall_err : array of float
        overall formal error in the VLBI solution
    num_obs : array of int
        number of obseravtions used in the solution
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(num_obs1, overall_err1, "r.", markersize=1, label="odd")
    ax.plot(num_obs2, overall_err2, "bx", markersize=1, label="even")

    N = np.arange(1, 1e6)
    sigma = 1. / np.sqrt(N)
    ax.plot(N, sigma, "k", label="$y=1.0/\\sqrt{N}$")

    ax.set_ylabel("Overall formal error (mas)")
    ax.set_yscale("log")
    ax.set_xlabel("$N_{observation}$")
    ax.set_xscale("log")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def comparison_plot(dat1, dat2, label1, label2, output):
    """Data comparison plot of two catalog

    Parameters
    ----------
    dat1/dat2 : array of float
        data from two catalogs
    label1/label2 : string
        label for x/y axis
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(dat1, dat2, "+", markersize=1)

    datmax = max(max(dat1), max(dat2))
    Nmax = 0.1
    while Nmax < datmax:
        Nmax *= 10
    N = np.arange(0, 1, 0.01) * Nmax
    ax.plot(N, N, "k", label="$y=N$")

    ax.set_xscale("log")
    ax.set_xlabel(label1)
    ax.set_yscale("log")
    ax.set_ylabel(label2)

    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def difference_vs_error_plot(diff_2, error_2, labelx, labely, output):
    """Data comparison of two catalog

    Parameters
    ----------
    diff_2/error_2 : array of float
        square of positional difference/combined error between two catalogs
    labelx/labely : string
        label for x/y axis
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(error_2, diff_2, "+", markersize=1)

    datmax = max(max(diff_2), max(error_2))
    Nmax = 0.1
    while Nmax < datmax:
        Nmax *= 10
    N = np.arange(0, 1, 0.01) * Nmax
    ax.plot(N, N, "k", label="$y=N$")

    ax.set_xscale("log")
    ax.set_xlabel(labelx)
    ax.set_yscale("log")
    ax.set_ylabel(labely)

    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def comparison_plot2(dat1_1, dat1_2, dat1_3,
                     dat2_1, dat2_2, dat2_3,
                     label1x, label1y, label2x, label2y, output):
    """Data comparison plot of two catalog

    Parameters
    ----------
    dat1/dat2 : array of float
        data from two catalogs
    label1/label2 : string
        label for x/y axis of figure 1 and 2
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.plot(dat1_1, dat1_2, "bo", markersize=2, label="WRMS")
    ax0.plot(dat1_1, dat1_2, "g", lw=1)
    ax0.plot(dat1_1, dat1_3, "r^", markersize=2, label="Inflated error")
    ax1.plot(dat2_1, dat2_2, "bo", markersize=2, label="WRMS")
    ax1.plot(dat2_1, dat2_2, "g", lw=1)
    ax1.plot(dat2_1, dat2_3, "r^", markersize=2, label="Inflated error")

    N = np.arange(min(min(dat1_1), min(dat2_1)), 1, 0.01)
    # ax0.plot(N, N, "k", label="$y=N$", lw=1)
    # ax1.plot(N, N, "k", label="$y=N$", lw=1)
    ax0.plot(N, N, "k", lw=0.5)
    ax1.plot(N, N, "k", lw=0.5)

    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_xlabel(label1x)
    # ax0.set_ylim([0, 2])
    ax0.set_ylabel(label1y)
    ax0.grid(True)
    ax0.legend()

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(label2x)
    # ax1.set_ylim([0, 2])
    ax1.set_ylabel(label2y)
    ax1.grid(True)
    ax1.legend()

    ax1.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def comparison_plot3(dat1_1, dat1_2, dat1_3, dat2_1, dat2_2, dat2_3,
                     label1x, label1y, label2x, label2y, output):
    """Data comparison plot of two catalog

    Parameters
    ----------
    dat1_1,2,3/dat2_1,2,3 : array of float
        data from two catalogs
    label1x,y/label2x,y : string
        label for x/y axis of figure 1 and 2
    output : string
        output filename with full path, ending with .eps or .png, etc.

    Returns
    ----------
    None
    """

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.plot(dat1_1, dat1_2, "bo", markersize=2, label="wrms noise")
    ax0.plot(dat1_1, dat1_3, "r^", markersize=2, label="median formal error")
    ax1.plot(dat2_1, dat2_2, "bo", markersize=2, label="wrms noise")
    ax1.plot(dat2_1, dat2_3, "r^", markersize=2, label="median formal error")

    # N = np.arange(0.01, 1, 0.01)
    # ax0.plot(N, N, "k", label="$y=N$", lw=1)
    # ax1.plot(N, N, "k", label="$y=N$", lw=1)

    ax0.set_xscale("log")
    ax0.set_xlabel(label1x)
    # ax0.set_yscale("log")
    ax0.set_ylabel(label1y)
    ax0.set_ylim([0, 0.25])
    ax0.grid(True)
    ax0.legend()

    ax1.set_xscale("log")
    # ax1.set_xlabel(label2x)
    # ax1.set_yscale("log")
    ax1.set_ylabel(label2y)
    ax1.set_ylim([0, 0.25])
    ax1.grid(True)
    ax1.legend()

    ax1.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def comparison_plot3_1(dat1, dat2, dat3, labelx, labely, output,
                       s=None, f=None, ylog=False):
    """Data comparison plot of two catalog

    Parameters
    ----------
    dat1,2,3 : array of float
        data from two catalogs
    label1x,y : string
        label for x/y axis
    output : string
        output filename with full path, ending with .eps or .png, etc.
    s : float
        scaling factor
    f : float
        noise floor
    ylog : boolean
        flag to determine if the scale of y-axis chose to be log

    Returns
    ----------
    None
    """

    fig, ax = plt.subplots()

    ax.plot(dat1, dat2, "bo", markersize=2, label="wrms noise")
    ax.plot(dat1, dat3, "r^", markersize=2, label="median formal error")

    if not s is None and not f is None:
        dat4 = np.sqrt((dat3 * s)**2 + f**2)
        ax.plot(dat1, dat4, "k", lw=0.5, label="Inflated error")

    ax.set_xscale("log")
    ax.set_xlabel(labelx)
    if ylog:
        ax.set_yscale("log")
    else:
        ax.set_ylim([0, 0.25])
    ax.set_ylabel(labely)
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()


def sf_nf_dec_plot(dec, sf_RA, nf_RA, sf_Dec, nf_Dec, output):
    """Data comparison plot of two catalog

    Parameters
    ----------

    Returns
    ----------
    None
    """

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    ax0.plot(dec, sf_RA, "bo-", markersize=5, label="R.A.")
    ax0.plot(dec, sf_Dec, "r^-", markersize=5, label="Dec.")
    ax1.plot(dec, nf_RA, "bo-", markersize=5, label="R.A.")
    ax1.plot(dec, nf_Dec, "r^-", markersize=5, label="Dec.")

    # ax0.set_xscale("log")
    ax0.set_xlabel("Declination (deg)")
    ax0.set_xlim([-90, 90])
    ax0.set_ylabel("Scale Factor")
    ax0.set_ylim([0.5, 2])
    ax0.grid(True)
    ax0.legend()

    # ax1.set_xscale("log")
    # ax1.set_xlabel(label2x)
    # ax1.set_yscale("log")
    ax1.set_ylabel("Noise Floor (mas)")
    ax1.set_ylim([0, 0.1])
    ax1.grid(True)
    ax1.legend()

    ax1.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    # plt.show()

    plt.savefig(output)
    plt.close()
# --------------------------------- END --------------------------------
