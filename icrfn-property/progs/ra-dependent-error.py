#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: dec-dependent-error.py
"""
Created on Sat May 18 10:59:56 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""


# -----------------------------  FUNCTIONS -----------------------------
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

# My modules
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3
from my_progs.catalog.read_gaia import read_dr2_iers
# from my_progs.catalog.smooth_err import smooth_by_dec


# Read the information from catalogs.

# ICRF1 catalog
icrf1 = read_icrf1()
icrf1.keep_columns(["ra", "ra_err", "dec_err", "pos_err"])
icrf1["ra_err"] = icrf1["ra_err"].to(u.uas)
icrf1["dec_err"] = icrf1["dec_err"].to(u.uas)
icrf1 = icrf1.filled()

# ICRF2 catalog
icrf2 = read_icrf2()
icrf2.keep_columns(["ra", "ra_err", "dec_err", "pos_err"])
icrf2["ra_err"] = icrf2["ra_err"].to(u.uas)
icrf2["dec_err"] = icrf2["dec_err"].to(u.uas)

# ICRF3 S/X catalog
icrf3sx = read_icrf3(wv="sx")
icrf3sx.keep_columns(["ra", "ra_err", "dec_err", "pos_err"])
icrf3sx["ra_err"] = icrf3sx["ra_err"].to(u.uas)
icrf3sx["dec_err"] = icrf3sx["dec_err"].to(u.uas)

# ICRF3 K catalog
icrf3k = read_icrf3(wv="k")
icrf3k.keep_columns(["ra", "ra_err", "dec_err", "pos_err"])
icrf3k["ra_err"] = icrf3k["ra_err"].to(u.uas)
icrf3k["dec_err"] = icrf3k["dec_err"].to(u.uas)

# ICRF3 X/Ka catalog
icrf3xka = read_icrf3(wv="xka")
icrf3xka.keep_columns(["ra", "ra_err", "dec_err", "pos_err"])
icrf3xka["ra_err"] = icrf3xka["ra_err"].to(u.uas)
icrf3xka["dec_err"] = icrf3xka["dec_err"].to(u.uas)


# Gaia DR2 aux_iers catalog
gdr2 = read_dr2_iers()
gdr2 = gdr2[gdr2["phot_g_mean_mag"] < 18.7].filled()
gdr2.keep_columns(["ra", "ra_err", "dec_err", "pos_err"])
gdr2["ra_err"] = gdr2["ra_err"].to(u.uas)
gdr2["dec_err"] = gdr2["dec_err"].to(u.uas)

# Bin the sources

binsize = 20

# ICRF1
icrf1_dec_bin = np.trunc(np.arange(len(icrf1)) / binsize)
icrf1.sort("ra")
icrf1_grouped = icrf1.group_by(icrf1_dec_bin)
icrf1_binned = icrf1_grouped.groups.aggregate(np.median)

binsize = 50

# ICRF2
icrf2_dec_bin = np.trunc(np.arange(len(icrf2)) / binsize)
icrf2.sort("ra")
icrf2_grouped = icrf2.group_by(icrf2_dec_bin)
icrf2_binned = icrf2_grouped.groups.aggregate(np.median)

# ICRF3 SX
icrf3sx_dec_bin = np.trunc(np.arange(len(icrf3sx)) / binsize)
icrf3sx.sort("ra")
icrf3sx_grouped = icrf3sx.group_by(icrf3sx_dec_bin)
icrf3sx_binned = icrf3sx_grouped.groups.aggregate(np.median)

binsize = 20

# ICRF3 K
icrf3k_dec_bin = np.trunc(np.arange(len(icrf3k)) / binsize)
icrf3k.sort("ra")
icrf3k_grouped = icrf3k.group_by(icrf3k_dec_bin)
icrf3k_binned = icrf3k_grouped.groups.aggregate(np.median)

# ICRF3 XKa
icrf3xka_dec_bin = np.trunc(np.arange(len(icrf3xka)) / binsize)
icrf3xka.sort("ra")
icrf3xka_grouped = icrf3xka.group_by(icrf3xka_dec_bin)
icrf3xka_binned = icrf3xka_grouped.groups.aggregate(np.median)

binsize = 50

# Gaia DR2
gdr2_dec_bin = np.trunc(np.arange(len(gdr2)) / binsize)
gdr2.sort("ra")
gdr2_grouped = gdr2.group_by(gdr2_dec_bin)
gdr2_binned = gdr2_grouped.groups.aggregate(np.median)


# Plot for median error
# fig, (ax0, ax1) = plt.subplots(figsize=(10, 4), ncols=2)
fig, (ax0, ax1) = plt.subplots(figsize=(4, 6), nrows=2, sharex=True)

# Right ascension
ax0.plot(gdr2_binned["ra"], gdr2_binned["ra_err"],
         "k-^", ms=2, label="Gaia DR2", lw=0.5)
ax0.plot(icrf1_binned["ra"], icrf1_binned["ra_err"],
         "g-v", ms=2, label="ICRF1", lw=0.5)
ax0.plot(icrf2_binned["ra"], icrf2_binned["ra_err"],
         "m-s", ms=2, label="ICRF2", lw=0.5)
ax0.plot(icrf3sx_binned["ra"], icrf3sx_binned["ra_err"],
         "r-x", ms=2, label="ICRF3 S/X", lw=0.5)
ax0.plot(icrf3k_binned["ra"], icrf3k_binned["ra_err"],
         "y-o", ms=2, label="ICRF3 K", lw=0.5)
ax0.plot(icrf3xka_binned["ra"], icrf3xka_binned["ra_err"],
         "b-*", ms=2, label="ICRF3 X/Ka", lw=0.5)

ax0.set_xlim([0, 360])
ax0.set_ylim([40, 6000])
ax0.set_xticks(np.arange(0, 361, 60))
# ax0.set_xlabel("Declination [$^{\\circ}$]", fontsize=10)
ax0.set_yscale("log")
# ax0.set_ylabel(
#     "Formal Error in $\\alpha^*$ [$\\mathrm{\\mu as}$]", fontsize=10)
ax0.set_ylabel(
    "$\\sigma_{\\alpha^*}$ ($\\mathrm{\\mu as}$)", fontsize=10)
# ax0.grid()
ax0.legend(fontsize="xx-small", loc="upper right")

ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")

ax0.xaxis.set_minor_locator(MultipleLocator(15))

# Declination
ax1.plot(gdr2_binned["ra"], gdr2_binned["dec_err"],
         "k-^", ms=2, label="Gaia DR2", lw=0.5)
ax1.plot(icrf1_binned["ra"], icrf1_binned["dec_err"],
         "g-v", ms=2, label="ICRF1", lw=0.5)
ax1.plot(icrf2_binned["ra"], icrf2_binned["dec_err"],
         "m-s", ms=2, label="ICRF2", lw=0.5)
ax1.plot(icrf3sx_binned["ra"], icrf3sx_binned["dec_err"],
         "r-x", ms=2, label="ICRF3 S/X", lw=0.5)
ax1.plot(icrf3k_binned["ra"], icrf3k_binned["dec_err"],
         "y-o", ms=2, label="ICRF3 K", lw=0.5)
ax1.plot(icrf3xka_binned["ra"], icrf3xka_binned["dec_err"],
         "b-*", ms=2, label="ICRF3 X/Ka", lw=0.5)

ax1.set_xlim([0, 360])
ax1.set_ylim([40, 6000])
ax1.set_xticks(np.arange(0, 361, 60))
ax1.set_xlabel("Right ascension ($^{\\circ}$)", fontsize=10)
ax1.set_yscale("log")
# ax1.set_ylabel("Formal Error in $\\delta$ [$\\mathrm{\\mu as}$]", fontsize=10)
ax1.set_ylabel("$\\sigma_\\delta$ ($\\mathrm{\\mu as}$)", fontsize=10)
# ax1.grid()
ax1.legend(fontsize="xx-small", loc="upper right")

ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")

ax1.xaxis.set_minor_locator(MultipleLocator(15))

plt.tight_layout()

# plt.savefig("../plots/pos-err-vs-dec.eps")
plt.savefig("../plots/pos-err-vs-ra_bgt.eps", hbox="tight")
# --------------------------------- END --------------------------------
