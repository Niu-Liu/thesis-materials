#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: g-limiting-magnitude-on-systematics.py
"""
Created on Fri Nov 15 14:19:27 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Check the global systematics of ICRF catalog as a function of g limiting magnitude
of the Gaia DR2 sample.
In fact, this is an attempt of checking the

"""

from astropy.table import setdiff
from astropy.table import Table, join, Column
from astropy import units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt

# My modules
from my_progs.catalog.vsh_deg2_cor import vsh_deg02_fitting, residual_calc02
from my_progs.catalog.pos_diff import nor_sep_calc, pos_diff_calc
from my_progs.catalog.pos_err import pos_err_calc
from my_progs.stat_func.rms_calc import rms_calc
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3
from my_progs.catalog.read_gaia import read_dr2_iers


# -----------------------------  FUNCTIONS -----------------------------
def outlier_elim(comsou):
    """Eliminate the outliers.
    """

    # Remove the outlier (consider the normalized separation)
    N = len(comsou)
    X0 = np.sqrt(np.log(N) * 2)

    mask = ((comsou["nor_sep"] <= X0)
            & (comsou["ang_sep"] < 10))

    # Table of a clean sample
    comcln = comsou[mask]
    N1 = len(comcln)
    outrate = (N - N1) / N

    return comcln, X0, N, N1, outrate


def vsh_estimate(comsou):
    """Estimate the VSH parameters from a sample.
    """

    # Transform columns into np.array
    dra = np.array(comsou["dra"])
    ddec = np.array(comsou["ddec"])
    dra_err = np.array(comsou["dra_err"])
    ddec_err = np.array(comsou["ddec_err"])
    ra_rad = np.array(comsou["ra_icrf"].to(u.radian))
    dec_rad = np.array(comsou["dec_icrf"].to(u.radian))
    dra_ddec_cov = np.array(comsou["dra_ddec_cov"])

    # Transformation parameters
    # l_max = 2
    w2_all, sig2_all, corrcoef2_all = vsh_deg02_fitting(
        dra, ddec, ra_rad, dec_rad, dra_err, ddec_err,
        cov=dra_ddec_cov, elim_flag="None")

    # mas -> uas
    w2 = w2_all * 1.e3
    sig2 = sig2_all * 1.e3

    return w2, sig2


# ----------------------------- MAIN -----------------------------
# Read Gaia DR2 IERS quasar data
gaiadr2 = read_dr2_iers()

# Read ICRF catalog
icrf1 = read_icrf1()
icrf2 = read_icrf2()
icrf3sx = read_icrf3(wv="sx")
icrf3k = read_icrf3(wv="k")
icrf3xka = read_icrf3(wv="xka")

# Crossmatch ICRF catalogs with Gaia DR2
comsou1 = join(icrf1, gaiadr2, keys="iers_name", table_names=["icrf", "gaia"])
comsou2 = join(icrf2, gaiadr2, keys="iers_name", table_names=["icrf", "gaia"])
comsou3sx = join(icrf3sx, gaiadr2, keys="iers_name", table_names=["icrf", "gaia"])
comsou3k = join(icrf3k, gaiadr2, keys="iers_name", table_names=["icrf", "gaia"])
comsou3xka = join(icrf3xka, gaiadr2, keys="iers_name", table_names=["icrf", "gaia"])

# Sort the table by G magnitude
comsou1.sort("phot_g_mean_mag")
comsou2.sort("phot_g_mean_mag")
comsou3sx.sort("phot_g_mean_mag")
comsou3k.sort("phot_g_mean_mag")
comsou3xka.sort("phot_g_mean_mag")

# Calculate the position offset
# ICRF1 - Gaia DR2
[dRA1, dDC1, dRAerr1, dDCerr1, dRAdDCcov1, angsep1, Xa1, Xd1, X1] = pos_diff_calc(
    comsou1["ra_icrf"], comsou1["ra_err_icrf"],
    comsou1["dec_icrf"], comsou1["dec_err_icrf"],
    comsou1["ra_dec_corr_icrf"], comsou1["ra_gaia"],
    comsou1["ra_err_gaia"], comsou1["dec_gaia"],
    comsou1["dec_err_gaia"], comsou1["ra_dec_corr_gaia"])

comsou1.add_columns([dRA1, dDC1, dRAerr1, dDCerr1, dRAdDCcov1, angsep1, Xa1, Xd1, X1],
                    names=["dra", "ddec", "dra_err", "ddec_err", "dra_ddec_cov",
                           "ang_sep", "nor_dra", "nor_ddec", "nor_sep"])

# ICRF2 - Gaia DR2
[dRA2, dDC2, dRAerr2, dDCerr2, dRAdDCcov2, angsep2, Xa2, Xd2, X2] = pos_diff_calc(
    comsou2["ra_icrf"], comsou2["ra_err_icrf"],
    comsou2["dec_icrf"], comsou2["dec_err_icrf"],
    comsou2["ra_dec_corr_icrf"], comsou2["ra_gaia"],
    comsou2["ra_err_gaia"], comsou2["dec_gaia"],
    comsou2["dec_err_gaia"], comsou2["ra_dec_corr_gaia"])

comsou2.add_columns([dRA2, dDC2, dRAerr2, dDCerr2, dRAdDCcov2, angsep2, Xa2, Xd2, X2],
                    names=["dra", "ddec", "dra_err", "ddec_err", "dra_ddec_cov",
                           "ang_sep", "nor_dra", "nor_ddec", "nor_sep"])

# ICRF3 SX - Gaia DR2
[dRA3sx, dDC3sx, dRAerr3sx, dDCerr3sx, dRAdDCcov3sx,
 angsep3sx, Xa3sx, Xd3sx, X3sx] = pos_diff_calc(
    comsou3sx["ra_icrf"], comsou3sx["ra_err_icrf"],
    comsou3sx["dec_icrf"], comsou3sx["dec_err_icrf"],
    comsou3sx["ra_dec_corr_icrf"], comsou3sx["ra_gaia"],
    comsou3sx["ra_err_gaia"], comsou3sx["dec_gaia"],
    comsou3sx["dec_err_gaia"], comsou3sx["ra_dec_corr_gaia"])

comsou3sx.add_columns([dRA3sx, dDC3sx, dRAerr3sx, dDCerr3sx, dRAdDCcov3sx,
                       angsep3sx, Xa3sx, Xd3sx, X3sx],
                      names=["dra", "ddec", "dra_err", "ddec_err", "dra_ddec_cov",
                             "ang_sep", "nor_dra", "nor_ddec", "nor_sep"])

# ICRF3 K - gaia DR2
[dRA3k, dDC3k, dRAerr3k, dDCerr3k, dRAdDCcov3k,
 angsep3k, Xa3k, Xd3k, X3k] = pos_diff_calc(
    comsou3k["ra_icrf"], comsou3k["ra_err_icrf"],
    comsou3k["dec_icrf"], comsou3k["dec_err_icrf"],
    comsou3k["ra_dec_corr_icrf"], comsou3k["ra_gaia"],
    comsou3k["ra_err_gaia"], comsou3k["dec_gaia"],
    comsou3k["dec_err_gaia"], comsou3k["ra_dec_corr_gaia"])

comsou3k.add_columns([dRA3k, dDC3k, dRAerr3k, dDCerr3k,
                      dRAdDCcov3k, angsep3k, Xa3k, Xd3k, X3k],
                     names=["dra", "ddec", "dra_err", "ddec_err", "dra_ddec_cov",
                            "ang_sep", "nor_dra", "nor_ddec", "nor_sep"])

# ICRF3 XKa - gaia DR2
[dRA3xka, dDC3xka, dRAerr3xka, dDCerr3xka, dRAdDCcov3xka,
 angsep3xka, Xa3xka, Xd3xka, X3xka] = pos_diff_calc(
    comsou3xka["ra_icrf"], comsou3xka["ra_err_icrf"],
    comsou3xka["dec_icrf"], comsou3xka["dec_err_icrf"],
    comsou3xka["ra_dec_corr_icrf"], comsou3xka["ra_gaia"],
    comsou3xka["ra_err_gaia"], comsou3xka["dec_gaia"],
    comsou3xka["dec_err_gaia"], comsou3xka["ra_dec_corr_gaia"])

comsou3xka.add_columns([dRA3xka, dDC3xka, dRAerr3xka, dDCerr3xka, dRAdDCcov3xka,
                        angsep3xka, Xa3xka, Xd3xka, X3xka],
                       names=["dra", "ddec", "dra_err", "ddec_err", "dra_ddec_cov",
                              "ang_sep", "nor_dra", "nor_ddec", "nor_sep"])

# Find the G magnitude of 100th source
icrf1g = comsou1[99]["phot_g_mean_mag"]
icrf2g = comsou2[99]["phot_g_mean_mag"]
icrf3sxg = comsou3sx[99]["phot_g_mean_mag"]
icrf3kg = comsou3k[99]["phot_g_mean_mag"]
icrf3xkag = comsou3xka[99]["phot_g_mean_mag"]
gmax0 = max([icrf1g, icrf2g, icrf3sxg, icrf3kg, icrf3xkag])

# Write this information in a text file
with open("../logs/g-limiting-mag-100th.log", "w") as flog:
    # Print to the screen
    print("Between ICRF1     and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf1g))
    print("Between ICRF2     and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf2g))
    print("Between ICRF3 SX  and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf3sxg))
    print("Between ICRF3 K   and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf3kg))
    print("Between ICRF3 XKa and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf3xkag))

    # Print into a log file
    print("Between ICRF1 and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf1g), file=flog)
    print("Between ICRF2 and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf2g), file=flog)
    print("Between ICRF3 SX and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf3sxg), file=flog)
    print("Between ICRF3 K and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf3kg), file=flog)
    print("Between ICRF3 XKa and Gaia DR2, the 100th source corresponds to g={:.2f} mag".format(
        icrf3xkag), file=flog)


# Create an array of maximum g magnitude
gmax0 = 17.5
gmaxbin = np.arange(gmax0, 21.01, 0.25)

# Write the output into several text file
# Statistics
# ICRF1
fsta1 = open("../logs/g-limiting-mag-stat-icrf1.log", "w")
print("# G-mag  NO.all  NO.cln  X0  rate", file=fsta1)
# ICRF2
fsta2 = open("../logs/g-limiting-mag-stat-icrf2.log", "w")
print("# G-mag  NO.all  NO.cln  X0  rate", file=fsta2)
# ICRF3 SX
fsta3sx = open("../logs/g-limiting-mag-stat-icrf3sx.log", "w")
print("# G-mag  NO.all  NO.cln  X0  rate", file=fsta3sx)
# ICRF3 K
fsta3k = open("../logs/g-limiting-mag-stat-icrf3k.log", "w")
print("# G-mag  NO.all  NO.cln  X0  rate", file=fsta3k)
# ICRF3 XKa
fsta3xka = open("../logs/g-limiting-mag-stat-icrf3xka.log", "w")
print("# G-mag  NO.all  NO.cln  X0  rate", file=fsta3xka)

# VSH parameters
# ICRF1
fvsh1 = open("../logs/g-limiting-mag-vsh02-icrf1.log", "w")
print("# G-mag  Glide(3)  Rotation(3)  Quadrupole(10)  +/-  \n"
      "# Unit: micro-arcsec", file=fvsh1)
# ICRF2
fvsh2 = open("../logs/g-limiting-mag-vsh02-icrf2.log", "w")
print("# G-mag  Glide(3)  Rotation(3)  Quadrupole(10)  +/-  \n"
      "# Unit: micro-arcsec", file=fvsh2)
# ICRF3 SX
fvsh3sx = open("../logs/g-limiting-mag-vsh02-icrf3sx.log", "w")
print("# G-mag  Glide(3)  Rotation(3)  Quadrupole(10)  +/-  \n"
      "# Unit: micro-arcsec", file=fvsh3sx)
# ICRF3 K
fvsh3k = open("../logs/g-limiting-mag-vsh02-icrf3k.log", "w")
print("# G-mag  Glide(3)  Rotation(3)  Quadrupole(10)  +/-  \n"
      "# Unit: micro-arcsec", file=fvsh3k)
# ICRF3 XKa
fvsh3xka = open("../logs/g-limiting-mag-vsh02-icrf3xka.log", "w")
print("# G-mag  Glide(3)  Rotation(3)  Quadrupole(10)  +/-  \n"
      "# Unit: micro-arcsec", file=fvsh3xka)

# Iterate in each bin to estimate the VSH parameters
for gmax in gmaxbin:
    # Get a sub-sample from common sources between ICRF1 and Gaia DR2
    mask1 = (comsou1["phot_g_mean_mag"] <= gmax)
    clnsou1, X0, N, N1, outrate = outlier_elim(comsou1[mask1])
    print("{:5.2f}  {:5d}  {:5d}  {:5.2f}  {:.2f}".format(gmax, N, N1, X0, outrate),
          file=fsta1)
    w2, sig2 = vsh_estimate(clnsou1)

    parline = "{:+6.0f}  " * 16
    errline = "{:6.0f}  " * 16
    print("{:5.2f}".format(gmax), parline.format(*w2), errline.format(*sig2), file=fvsh1)

    # Same trick between ICRF2 and Gaia DR2
    mask2 = (comsou2["phot_g_mean_mag"] <= gmax)
    clnsou2, X0, N, N1, outrate = outlier_elim(comsou2[mask2])
    print("{:5.2f}  {:5d}  {:5d}  {:5.2f}  {:.2f}".format(gmax, N, N1, X0, outrate),
          file=fsta2)
    w2, sig2 = vsh_estimate(clnsou2)
    print("{:5.2f}".format(gmax), parline.format(*w2), errline.format(*sig2), file=fvsh2)

    # Same trick between ICRF3 SX and Gaia DR2
    mask3sx = (comsou3sx["phot_g_mean_mag"] <= gmax)
    clnsou3sx, X0, N, N1, outrate = outlier_elim(comsou3sx[mask3sx])
    print("{:5.2f}  {:5d}  {:5d}  {:5.2f}  {:.2f}".format(gmax, N, N1, X0, outrate),
          file=fsta3sx)
    w2, sig2 = vsh_estimate(clnsou3sx)
    print("{:5.2f}".format(gmax), parline.format(*w2), errline.format(*sig2), file=fvsh3sx)

    # Same trick between ICRF3 K and Gaia DR2
    mask3k = (comsou3k["phot_g_mean_mag"] <= gmax)
    clnsou3k, X0, N, N1, outrate = outlier_elim(comsou3k[mask3k])
    print("{:5.2f}  {:5d}  {:5d}  {:5.2f}  {:.2f}".format(gmax, N, N1, X0, outrate),
          file=fsta3k)
    w2, sig2 = vsh_estimate(clnsou3k)
    print("{:5.2f}".format(gmax), parline.format(*w2), errline.format(*sig2), file=fvsh3k)

    # Same trick between ICRF3 XKa and Gaia DR2
    mask3xka = (comsou3xka["phot_g_mean_mag"] <= gmax)
    clnsou3xka, X0, N, N1, outrate = outlier_elim(comsou3xka[mask3xka])
    print("{:5.2f}  {:5d}  {:5d}  {:5.2f}  {:.2f}".format(gmax, N, N1, X0, outrate),
          file=fsta3xka)
    w2, sig2 = vsh_estimate(clnsou3xka)
    print("{:5.2f}".format(gmax), parline.format(*w2), errline.format(*sig2), file=fvsh3xka)

# Close all the open files
fsta1.close()
fsta2.close()
fsta3sx.close()
fsta3k.close()
fsta3xka.close()
fvsh1.close()
fvsh2.close()
fvsh3sx.close()
fvsh3k.close()
fvsh3xka.close()
# --------------- END - -------------------------------
