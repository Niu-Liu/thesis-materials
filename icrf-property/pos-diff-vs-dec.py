#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: pos-diff-vs-dec.py
"""
Created on Sat May 18 10:33:38 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

# My modules
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3
from my_progs.catalog.read_gaia import read_dr2_iers
from my_progs.catalog.pos_diff import radio_cat_diff_calc
from my_progs.catalog.smoothed_error import smooth_by_dec


# Read the catalogs.
# ICRF1 catalog
icrf1 = read_icrf1()

# ICRF2 catalog
icrf2 = read_icrf2()

# ICRF3 S/X catalog
icrf3sx = read_icrf3(wv="sx")

# ICRF3 K catalog
icrf3k = read_icrf3(wv="k")

# ICRF3 X/Ka catalog
icrf3xka = read_icrf3(wv="xka")

# Gaia DR2 aux_iers catalog
gdr2 = read_dr2_iers()
# Consider only quasars with G<18.7
mask = (gdr2["phot_g_mean_mag"] < 18.7)
gdr2 = gdr2[mask].filled()

# Calculate the positional offsets.

# ICRF1 - Gaia DR2
icrf1_gdr2 = radio_cat_diff_calc(icrf1, gdr2, "iers_name")
icrf1_gdr2.keep_columns(["dec", "dra", "ddec"])
icrf1_gdr2["dra"] = icrf1_gdr2["dra"].to(u.uas)
icrf1_gdr2["ddec"] = icrf1_gdr2["ddec"].to(u.uas)
icrf1_gdr2.filled()

# ICRF2 - Gaia DR2
icrf2_gdr2 = radio_cat_diff_calc(icrf2, gdr2, "iers_name")
icrf2_gdr2.keep_columns(["dec", "dra", "ddec"])
icrf2_gdr2["dra"] = icrf2_gdr2["dra"].to(u.uas)
icrf2_gdr2["ddec"] = icrf2_gdr2["ddec"].to(u.uas)
icrf2_gdr2.filled()

# ICRF3 SX - Gaia DR2
icrf3sx_gdr2 = radio_cat_diff_calc(icrf3sx, gdr2, "iers_name")
icrf3sx_gdr2.keep_columns(["dec", "dra", "ddec"])
icrf3sx_gdr2["dra"] = icrf3sx_gdr2["dra"].to(u.uas)
icrf3sx_gdr2["ddec"] = icrf3sx_gdr2["ddec"].to(u.uas)
icrf3sx_gdr2.filled()

# ICRF3 K - Gaia DR2
icrf3k_gdr2 = radio_cat_diff_calc(icrf3k, gdr2, "iers_name")
icrf3k_gdr2.keep_columns(["dec", "dra", "ddec"])
icrf3k_gdr2["dra"] = icrf3k_gdr2["dra"].to(u.uas)
icrf3k_gdr2["ddec"] = icrf3k_gdr2["ddec"].to(u.uas)
icrf3k_gdr2.filled()

# ICRF3 Xka Gaia DR2
icrf3xka_gdr2 = radio_cat_diff_calc(icrf3xka, gdr2, "iers_name")
icrf3xka_gdr2.keep_columns(["dec", "dra", "ddec"])
icrf3xka_gdr2["dra"] = icrf3xka_gdr2["dra"].to(u.uas)
icrf3xka_gdr2["ddec"] = icrf3xka_gdr2["ddec"].to(u.uas)
icrf3xka_gdr2.filled()


# First, I try with binning sources with a size of 50 sources.
# binsize1 = 20
# binsize2 = 50
binsize1 = 10
binsize2 = 25

# ICRF1
icrf1_dec_bin = np.trunc(np.arange(len(icrf1_gdr2)) / binsize1)
icrf1_gdr2.sort("dec")
icrf1_grouped = icrf1_gdr2.group_by(icrf1_dec_bin)
icrf1_binned = icrf1_grouped.groups.aggregate(np.median)

# ICRF2
icrf2_dec_bin = np.trunc(np.arange(len(icrf2_gdr2)) / binsize2)
icrf2_gdr2.sort("dec")
icrf2_grouped = icrf2_gdr2.group_by(icrf2_dec_bin)
icrf2_binned = icrf2_grouped.groups.aggregate(np.median)

# ICRF3 SX
icrf3sx_dec_bin = np.trunc(np.arange(len(icrf3sx_gdr2)) / binsize2)
icrf3sx_gdr2.sort("dec")
icrf3sx_grouped = icrf3sx_gdr2.group_by(icrf3sx_dec_bin)
icrf3sx_binned = icrf3sx_grouped.groups.aggregate(np.median)

# ICRF3 K
icrf3k_dec_bin = np.trunc(np.arange(len(icrf3k_gdr2)) / binsize1)
icrf3k_gdr2.sort("dec")
icrf3k_grouped = icrf3k_gdr2.group_by(icrf3k_dec_bin)
icrf3k_binned = icrf3k_grouped.groups.aggregate(np.median)

# ICRF3 XKa
icrf3xka_dec_bin = np.trunc(np.arange(len(icrf3xka_gdr2)) / binsize1)
icrf3xka_gdr2.sort("dec")
icrf3xka_grouped = icrf3xka_gdr2.group_by(icrf3xka_dec_bin)
icrf3xka_binned = icrf3xka_grouped.groups.aggregate(np.median)

# Plot

minorLocator1 = MultipleLocator(100)
minorLocator2 = MultipleLocator(15)

# Plot for median error
# fig, (ax0, ax1) = plt.subplots(figsize=(10, 4), ncols=2)
# fig, (ax0, ax1) = plt.subplots(figsize=(10, 4), nrows=2)
fig, (ax0, ax1) = plt.subplots(figsize=(10, 4), ncols=2, sharex=True)

ax0.plot(icrf1_binned["dec"], icrf1_binned["dra"], "g-v", ms=2, label="ICRF1", lw=0.5)

ax0.plot(icrf2_binned["dec"], icrf2_binned["dra"], "m-s", ms=2, label="ICRF2", lw=0.5)

ax0.plot(icrf3sx_binned["dec"], icrf3sx_binned["dra"],
         "r-x", ms=2, label="ICRF3 S/X", lw=0.5)

ax0.plot(icrf3k_binned["dec"], icrf3k_binned["dra"],
         "y-o", ms=2, label="ICRF3 K", lw=0.5)

ax0.plot(icrf3xka_binned["dec"], icrf3xka_binned["dra"],
         "b-*", ms=2, label="ICRF3 X/Ka", lw=0.5)

# ax0.set_ylim([-700, 700])
ax0.set_ylim([-700, 700])
ax0.set_xlim([-90, 90])
ax0.set_xticks(np.arange(-90, 91, 30))

ax0.set_xlabel("Declination ($^{\\circ}$)", fontsize=10)
# ax0.set_ylabel("Difference in $\\alpha^*$ [$\\mathrm{\\mu as}$]", fontsize=10)
ax0.set_ylabel(
    "$\\Delta\\alpha^*$ ($\\mathrm{\\mu as}$)", fontsize=10)
ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")

ax0.yaxis.set_minor_locator(minorLocator1)
ax0.xaxis.set_minor_locator(minorLocator2)

# ax0.grid()
ax0.legend(fontsize="xx-small")

ax1.plot(icrf1_binned["dec"], icrf1_binned["ddec"], "g-v", ms=2, label="ICRF1", lw=0.5)

ax1.plot(icrf2_binned["dec"], icrf2_binned["ddec"], "m-s", ms=2, label="ICRF2", lw=0.5)

ax1.plot(icrf3sx_binned["dec"], icrf3sx_binned["ddec"],
         "r-x", ms=2, label="ICRF3 S/X", lw=0.5)

ax1.plot(icrf3k_binned["dec"], icrf3k_binned["ddec"],
         "y-o", ms=2, label="ICRF3 K", lw=0.5)

ax1.plot(icrf3xka_binned["dec"], icrf3xka_binned["ddec"],
         "b-*", ms=2, label="ICRF3 X/Ka", lw=0.5)

ax1.set_ylim([-700, 700])
ax1.set_xlim([-90, 90])
ax1.set_xticks(np.arange(-90, 91, 30))

ax1.set_xlabel("Declination ($^{\\circ}$)", fontsize=10)
# ax1.set_ylabel("Difference in $\\delta$ [$\\mathrm{\\mu as}$]", fontsize=10)
ax1.set_ylabel("$\\Delta\\delta$ ($\\mathrm{\\mu as}$)", fontsize=10)
ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")

ax1.yaxis.set_minor_locator(minorLocator1)
ax1.xaxis.set_minor_locator(minorLocator2)

# ax1.grid()
ax1.legend(fontsize="xx-small")

plt.tight_layout()
# plt.savefig("../plots/pos-diff-vs-dec.eps")
plt.savefig("../plots/pos-diff-vs-dec_bgt.eps", hbox="tight")
# plt.savefig("../plots/fig04.eps")
# --------------------------------- END --------------------------------
