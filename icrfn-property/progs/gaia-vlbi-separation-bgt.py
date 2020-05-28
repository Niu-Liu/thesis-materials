#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: gaia-vlbi-separation.py
"""
Created on Tue May 14 16:41:19 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""


from astropy.table import Table, join, Column
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

# My modules
from my_progs.catalog.vsh_deg2_cor import vsh_deg02_fitting, residual_calc02
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3
from my_progs.catalog.read_gaia import read_dr2_iers
from my_progs.catalog.pos_diff import radio_cat_diff_calc

# -----------------------------  MAIN -----------------------------
# ======================================================================
# Read the catalog
# ======================================================================
# Read Gaia DR2 IERS quasar data
gdr2 = read_dr2_iers()
gdr2 = gdr2[gdr2["phot_g_mean_mag"] < 18.7].filled()

# Read ICRF1 catalog
icrf1 = read_icrf1()
icrf1_type = Table(icrf1)
icrf1_type.keep_columns(["iers_name", "type"])

# Read ICRF2 catalog
icrf2 = read_icrf2()
icrf2_type = Table(icrf2)
icrf2_type.keep_columns(["iers_name", "type"])

# Read ICRF3 S/X catalog
icrf3sx = read_icrf3(wv="sx")
icrf3sx_type = Table(icrf3sx)
icrf3sx_type.keep_columns(["iers_name", "type"])

# Read ICRF3 K catalog
icrf3k = read_icrf3(wv="k")
icrf3k_type = Table(icrf3k)
icrf3k_type.keep_columns(["iers_name", "type"])

# Read ICRF3 S/X catalog
icrf3xka = read_icrf3(wv="xka")
icrf3xka_type = Table(icrf3xka)
icrf3xka_type.keep_columns(["iers_name", "type"])

# ======================================================================
# Calculate the position difference
# ======================================================================
# ICRF1 - Gaia DR2
icrf1_gdr2 = radio_cat_diff_calc(icrf1, gdr2, "iers_name", ["icrf1", "gdr2"])
icrf1_gdr2 = join(icrf1_gdr2, icrf1_type, keys="iers_name")

# ICRF2 - Gaia DR2
icrf2_gdr2 = radio_cat_diff_calc(icrf2, gdr2, "iers_name", ["icrf2", "gdr2"])
icrf2_gdr2 = join(icrf2_gdr2, icrf2_type, keys="iers_name")

# ICRF3SX - Gaia DR2
icrf3sx_gdr2 = radio_cat_diff_calc(
    icrf3sx, gdr2, "iers_name", ["icrf3sx", "gdr2"])
icrf3sx_gdr2 = join(icrf3sx_gdr2, icrf3sx_type, keys="iers_name")

# ICRF3K - Gaia DR2
icrf3k_gdr2 = radio_cat_diff_calc(
    icrf3k, gdr2, "iers_name", ["icrf3k", "gdr2"])
icrf3k_gdr2 = join(icrf3k_gdr2, icrf3k_type, keys="iers_name")

# ICRF3XKa - Gaia DR2
icrf3xka_gdr2 = radio_cat_diff_calc(
    icrf3xka, gdr2, "iers_name", ["icrf3xka", "gdr2"])
icrf3xka_gdr2 = join(icrf3xka_gdr2, icrf3xka_type, keys="iers_name")

# Unit transformation
icrf1_gdr2["ang_sep"] = icrf1_gdr2["ang_sep"].to(u.uas)
icrf2_gdr2["ang_sep"] = icrf2_gdr2["ang_sep"].to(u.uas)
icrf3sx_gdr2["ang_sep"] = icrf3sx_gdr2["ang_sep"].to(u.uas)
icrf3k_gdr2["ang_sep"] = icrf3k_gdr2["ang_sep"].to(u.uas)
icrf3xka_gdr2["ang_sep"] = icrf3xka_gdr2["ang_sep"].to(u.uas)

# Calculate X0
N1 = len(icrf1_gdr2)
N2 = len(icrf2_gdr2)
N3sx = len(icrf3sx_gdr2)
N3k = len(icrf3k_gdr2)
N3xka = len(icrf3xka_gdr2)

X0_1 = np.sqrt(2 * np.log(N1))
X0_2 = np.sqrt(2 * np.log(N2))
X0_3sx = np.sqrt(2 * np.log(N3sx))
X0_3k = np.sqrt(2 * np.log(N3k))
X0_3xka = np.sqrt(2 * np.log(N3xka))

print("X0 for each pair is ")
print("ICRF1     - Gaia DR2:  {:.2f}".format(X0_1))
print("ICRF2     - Gaia DR2:  {:.2f}".format(X0_2))
print("ICRF3 SX  - Gaia DR2:  {:.2f}".format(X0_3sx))
print("ICRF3 K   - Gaia DR2:  {:.2f}".format(X0_3k))
print("ICRF3 XKa - Gaia DR2:  {:.2f}".format(X0_3xka))

# ======================================================================
# Separate sources into groups of different categories
# ======================================================================
# Seperate the ICRF1 sources into the "defining", "Candidate", "Other".
# Defining sources
mask_def = (icrf1_gdr2["type"] == "D")
com_def_icrf1 = icrf1_gdr2[mask_def]

# Candiate sources
mask_can = (icrf1_gdr2["type"] == "C")
com_can_icrf1 = icrf1_gdr2[mask_can]

# Other sources
mask_oth = (icrf1_gdr2["type"] == "O")
com_oth_icrf1 = icrf1_gdr2[mask_oth]

# Seperate the ICRF2 sources into the "defining", "VCS-only" and "Non-VCS".
# Defining sources
mask_def = (icrf2_gdr2["type"] == "D")
com_def_icrf2 = icrf2_gdr2[mask_def]

# VCS sources
mask_vcs = (icrf2_gdr2["type"] == "V")
com_vcs_icrf2 = icrf2_gdr2[mask_vcs]

# Non-VCS sources
mask_nvc = (icrf2_gdr2["type"] == "N")
com_nvc_icrf2 = icrf2_gdr2[mask_nvc]

# Seperate the  ICRF3 S/X sources into the "defining" and "other".
# Defining sources
mask_def = (icrf3sx_gdr2["type"] == "D")
com_def_icrf3sx = icrf3sx_gdr2[mask_def]

# Other sources
mask_oth = (icrf3sx_gdr2["type"] != "D")
com_oth_icrf3sx = icrf3sx_gdr2[mask_oth]

# Seperate the  ICRF3 K sources into the "defining" and "other".
# Defining sources
mask_def = (icrf3k_gdr2["type"] == "D")
com_def_icrf3k = icrf3k_gdr2[mask_def]

# Other sources
mask_oth = (icrf3k_gdr2["type"] != "D")
com_oth_icrf3k = icrf3k_gdr2[mask_oth]

# Seperate the  ICRF3 X/Ka sources into the "defining" and "other".
# Defining sources
mask_def = (icrf3xka_gdr2["type"] == "D")
com_def_icrf3xka = icrf3xka_gdr2[mask_def]

# Other sources
mask_oth = (icrf3xka_gdr2["type"] != "D")
com_oth_icrf3xka = icrf3xka_gdr2[mask_oth]

# ======================================================================
# Make the plot
# ======================================================================
# Angular separation vs. normalized separation
fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(figsize=(6, 10), nrows=5,
                                              sharex=True)

# ICRF1
ax0.plot(com_can_icrf1["nor_sep"], com_can_icrf1["ang_sep"],
         "b*", ms=1, label="Candidate")
ax0.plot(com_oth_icrf1["nor_sep"], com_oth_icrf1["ang_sep"],
         "g^", ms=1, label="Other")
ax0.plot(com_def_icrf1["nor_sep"], com_def_icrf1["ang_sep"],
         "rx", ms=1, label="Defining")

ax0.hlines(10000, 0.01, 400, colors="k", linestyles="dashed", lw=0.1)
ax0.vlines(X0_1, 10, 600000, colors="k", linestyles="dashed", lw=0.1)

ax0.set_xlim([0.01, 400])
ax0.set_ylim([10, 600000])

ax0.set_xscale("log")
ax0.set_yscale("log")

# ax0.set_ylabel("Normalized separation $X$")
# ax0.set_xlabel("$X$")
ax0.set_ylabel("$\\rho$ ($\\mathrm{\\mu as}$)", fontsize=15)

ax0.text(30, 100, "ICRF1", fontsize=15)

ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")

ax0.legend(loc="upper left", fontsize="x-small")

# ICRF2
ax1.plot(com_vcs_icrf2["nor_sep"], com_vcs_icrf2["ang_sep"],
         "b*", ms=1, label="VCS-only")
ax1.plot(com_nvc_icrf2["nor_sep"], com_nvc_icrf2["ang_sep"],
         "g^", ms=1, label="Non-VCS")
ax1.plot(com_def_icrf2["nor_sep"], com_def_icrf2["ang_sep"],
         "rx", ms=1, label="Defining")

ax1.hlines(10000, 0.01, 400, colors="k", linestyles="dashed", lw=0.1)
ax1.vlines(X0_2, 10, 600000, colors="k", linestyles="dashed", lw=0.1)

ax1.set_xscale("log")
ax1.set_yscale("log")

ax1.set_ylabel("$\\rho$ ($\\mathrm{\\mu as}$)", fontsize=15)

ax1.text(30, 100, "ICRF2", fontsize=15)

# ax1.set_xlim([0.02, 150])
ax1.set_ylim([10, 600000])
ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")

ax1.legend(loc="upper left", fontsize="x-small")

# ICRF3 SX
ax2.plot(com_oth_icrf3sx["nor_sep"],
         com_oth_icrf3sx["ang_sep"], "b*", ms=1, label="Other")
ax2.plot(com_def_icrf3sx["nor_sep"],
         com_def_icrf3sx["ang_sep"], "rx", ms=1, label="Defining")
ax2.hlines(10000, 0.01, 400, colors="k", linestyles="dashed", lw=0.1)
ax2.vlines(X0_3sx, 10, 600000, colors="k", linestyles="dashed", lw=0.1)

ax2.set_ylabel("$\\rho$ ($\\mathrm{\\mu as}$)", fontsize=15)

ax2.text(30, 100, "ICRF3 S/X", fontsize=15)

ax2.set_xscale("log")
ax2.set_yscale("log")

# ax2.set_xlim([0.02, 600])
ax2.set_ylim([10, 600000])

ax2.xaxis.set_ticks_position("both")
ax2.yaxis.set_ticks_position("both")

ax2.legend(loc="upper left", fontsize="x-small")

# ICRF3 K
ax3.plot(com_oth_icrf3k["nor_sep"],
         com_oth_icrf3k["ang_sep"], "b*", ms=1, label="Other")
ax3.plot(com_def_icrf3k["nor_sep"],
         com_def_icrf3k["ang_sep"], "rx", ms=1, label="Defining")
ax3.hlines(10000, 0.01, 400, colors="k", linestyles="dashed", lw=0.1)
ax3.vlines(X0_3k, 10, 600000, colors="k", linestyles="dashed", lw=0.1)

ax3.set_xscale("log")
ax3.set_yscale("log")

ax3.set_ylabel("$\\rho$ ($\\mathrm{\\mu as}$)", fontsize=15)

ax3.text(30, 100, "ICRF3 K", fontsize=15)

# ax3.set_xlim([0.08, 100])
ax3.set_ylim([10, 600000])

ax3.xaxis.set_ticks_position("both")
ax3.yaxis.set_ticks_position("both")

ax3.legend(loc="upper left", fontsize="x-small")

# ICRF3 X/Ka
ax4.plot(com_oth_icrf3xka["nor_sep"], com_oth_icrf3xka["ang_sep"],
         "b*", ms=1, label="Other")
ax4.plot(com_def_icrf3xka["nor_sep"], com_def_icrf3xka["ang_sep"],
         "rx", ms=1, label="Defining")
ax4.hlines(10000, 0.01, 400, colors="k", linestyles="dashed", lw=0.1)
ax4.vlines(X0_3xka, 10, 600000, colors="k", linestyles="dashed", lw=0.1)

ax4.set_xscale("log")
ax4.set_yscale("log")

ax4.set_xlabel("Normalized seperation $X$", fontsize=15)
ax4.set_ylabel("$\\rho$ ($\\mathrm{\\mu as}$)", fontsize=15)

ax4.text(30, 100, "ICRF3 X/Ka", fontsize=15)

# ax4.set_xlim([0.08, 100])
ax4.set_ylim([10, 600000])

ax4.xaxis.set_ticks_position("both")
ax4.yaxis.set_ticks_position("both")

ax4.legend(loc="upper left", fontsize="x-small")

plt.subplots_adjust()
plt.tight_layout()
plt.savefig("../plots/gaia-vlbi-separation_bgt.eps", hbox="tight")
# plt.savefig("../plots/fig05.eps")
# --------------------------------- END --------------------------------
