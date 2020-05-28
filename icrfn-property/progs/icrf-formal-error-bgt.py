#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: icrf-formal-error.py
"""
Created on Sat May 18 11:30:42 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy import units as u
from my_progs.catalog.pos_err import pos_err_calc
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
# My modules
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3


# -----------------------------  FUNCTIONS -----------------------------
# I looked at the formal error of the ICRF catalog and Gaia DR2 to see if the formal errors are Gaussian-like.

# Read these ICRF catalogs from data files.

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

# Gaia DR2
dr2_file = ("/Users/Neo/Astronomy/data/catalogs/"
            "gaia/dr2/gaiadr2_iers.fits")


gdr2 = Table.read(dr2_file)
gadr2 = gdr2[gdr2["phot_g_mean_mag"] < 18.7]

pos_err = pos_err_calc(
    gdr2["ra_error"], gdr2["dec_error"], gdr2["ra_dec_corr"])

gdr2.add_column(pos_err, name="pos_err", index=6)

# Unit transformation
icrf1["pos_err"] = icrf1["pos_err"] * 1.e3
icrf2["pos_err"] = icrf2["pos_err"] * 1.e3
icrf3sx["pos_err"] = icrf3sx["pos_err"] * 1.e3
icrf3k["pos_err"] = icrf3k["pos_err"] * 1.e3
icrf3xka["pos_err"] = icrf3xka["pos_err"] * 1.e3
gdr2["pos_err"] = gdr2["pos_err"] * 1.e3


# Then I plot the formal against the number of observed delay and number of observed sessions.

# For ICRF1 catalog, there are three classifications of sources:
#
# 1. Defining
# 2. Candidates
# 3. Other
#
mask1 = icrf1["type"] == "D"
mask2 = icrf1["type"] == "C"
mask3 = icrf1["type"] == "O"

icrf1_def = icrf1[mask1]
icrf1_can = icrf1[mask2]
icrf1_oth = icrf1[mask3]

# For ICRF2 catalog, there are also three classifications of sources:
#
# 1. Defining
# 2. VCS-only
# 3. Non-VCS
#
mask1 = icrf2["type"] == "D"
mask2 = icrf2["type"] == "V"
mask3 = icrf2["type"] == "N"

icrf2_def = icrf2[mask1]
icrf2_vcs = icrf2[mask2]
icrf2_non = icrf2[mask3]

# For ICRF3 catalog, there are only two classifications of sources:
#
# 1. Defining
# 2. Other

mask1 = icrf3sx["type"] == "D"
mask2 = icrf3sx["type"] != "D"

icrf3sx_def = icrf3sx[mask1]
icrf3sx_oth = icrf3sx[mask2]

mask1 = icrf3k["type"] == "D"
mask2 = icrf3k["type"] != "D"

icrf3k_def = icrf3k[mask1]
icrf3k_oth = icrf3k[mask2]

mask1 = icrf3xka["type"] == "D"
mask2 = icrf3xka["type"] != "D"

icrf3xka_def = icrf3xka[mask1]
icrf3xka_oth = icrf3xka[mask2]

# Plot

fig, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(
    figsize=(10, 12), ncols=2, nrows=3)

# Gaia DR2
ax0.plot(gdr2["astrometric_n_obs_al"], gdr2["pos_err"], "bx", ms=1)

ax0.set_yscale("log")
ax0.set_xscale("log")
ax0.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ Gaia DR2 ($\\mathrm{\\mu as}$)", fontsize=12)
ax0.set_xlabel("astrometric_n_obs_al", fontsize=12)

# Gaussian-like formal
y = 1 / np.sqrt(gdr2["astrometric_n_obs_al"]) * 1e3
ax0.plot(gdr2["astrometric_n_obs_al"], y, "k-", label="", lw=0.1)

ax0.set_ylim([10, 10000])
ax0.yaxis.set_ticks_position("both")
ax0.xaxis.set_ticks_position("both")

# ICRF1
ax1.plot(icrf1_def["nb_del"], icrf1_def["pos_err"], "rx", label="Defining", ms=1)
ax1.plot(icrf1_can["nb_del"], icrf1_can["pos_err"], "b*", label="Candidate", ms=1)
ax1.plot(icrf1_oth["nb_del"], icrf1_oth["pos_err"], "g^", label="Other", ms=1)

# Gaussian-like formal
y = 3 / np.sqrt(icrf1["nb_del"]) * 1e3
ax1.plot(icrf1["nb_del"], y, "k-", label="", lw=0.1)

ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_ylim([100, 2000000])
ax1.yaxis.set_ticks_position("both")
ax1.xaxis.set_ticks_position("both")
ax1.legend(fontsize="x-small")
# ax1.text(100, 800, "ICRF1")
ax1.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF1 ($\\mathrm{\\mu as}$)", fontsize=12)
ax1.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF2
ax2.plot(icrf2_vcs["nb_del"], icrf2_vcs["pos_err"], "b*", label="VCS", ms=1)
ax2.plot(icrf2_non["nb_del"], icrf2_non["pos_err"], "g^", label="Non-VCS", ms=1)
ax2.plot(icrf2_def["nb_del"], icrf2_def["pos_err"], "rx", label="Defining", ms=1)

# Gaussian-like formal
y = 2 / np.sqrt(icrf2["nb_del"]) * 1e3
ax2.plot(icrf2["nb_del"], y, "k-", label="", lw=0.1)

ax2.set_yscale("log")
ax2.set_xscale("log")
ax2.set_ylim([20, 500000])
ax2.yaxis.set_ticks_position("both")
ax2.xaxis.set_ticks_position("both")
ax2.legend(fontsize="x-small")

ax2.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF2 ($\\mathrm{\\mu as}$)", fontsize=12)
ax2.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF3 SX
ax3.plot(icrf3sx_oth["nb_del"], icrf3sx_oth["pos_err"], "b*", label="Other", ms=1)
ax3.plot(icrf3sx_def["nb_del"], icrf3sx_def["pos_err"], "rx", label="Defining", ms=1)

# Gaussian-like formal
y = 2 / np.sqrt(icrf3sx["nb_del"]) * 1e3
ax3.plot(icrf3sx["nb_del"], y, "k-", label="", lw=0.1)

ax3.set_ylim([20, 500000])

ax3.set_yscale("log")
ax3.set_xscale("log")
ax3.yaxis.set_ticks_position("both")
ax3.xaxis.set_ticks_position("both")
ax3.legend(fontsize="x-small")

ax3.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF3 S/X ($\\mathrm{\\mu as}$)", fontsize=12)
ax3.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF3 K
ax4.plot(icrf3k_def["nb_del"], icrf3k_def["pos_err"], "rx", label="Defining", ms=1)
ax4.plot(icrf3k_oth["nb_del"], icrf3k_oth["pos_err"], "b*", label="Other", ms=1)

# Gaussian-like formal
y = 2 / np.sqrt(icrf3k["nb_del"]) * 1e3
ax4.plot(icrf3k["nb_del"], y, "k-", label="", lw=0.1)

ax4.set_ylim([20, 30000])
ax4.set_yscale("log")
ax4.set_xscale("log")
ax4.yaxis.set_ticks_position("both")
ax4.xaxis.set_ticks_position("both")
ax4.legend(fontsize="x-small")

ax4.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF3 K ($\\mathrm{\\mu as}$)", fontsize=12)
ax4.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF3 XKa
ax5.plot(icrf3xka_def["nb_del"],
         icrf3xka_def["pos_err"], "rx", label="Defining", ms=1)
ax5.plot(icrf3xka_oth["nb_del"], icrf3xka_oth["pos_err"], "b*", label="Other", ms=1)


# Gaussian-like formal
y = 1 / np.sqrt(icrf3xka["nb_del"]) * 1e3
ax5.plot(icrf3xka["nb_del"], y, "k-", label="", lw=0.1)

ax5.set_ylim([30, 200000])
ax5.set_yscale("log")
ax5.set_xscale("log")
ax5.yaxis.set_ticks_position("both")
ax5.xaxis.set_ticks_position("both")
ax5.legend(fontsize="x-small")

ax5.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF3 X/Ka ($\\mathrm{\\mu as}$)", fontsize=12)
ax5.set_xlabel("$N_{\\rm delay}$", fontsize=12)

plt.tight_layout()
plt.savefig("../plots/error-vs-nobs_bgt.eps", hbox="tight")
plt.close()

# ====================  South-North ============================
icrf1_nor = icrf1[icrf1["dec"] >= 0]
icrf1_sou = icrf1[icrf1["dec"] < 0]

icrf2_nor = icrf2[icrf2["dec"] >= 0]
icrf2_sou = icrf2[icrf2["dec"] < 0]

icrf3sx_nor = icrf3sx[icrf3sx["dec"] >= 0]
icrf3sx_sou = icrf3sx[icrf3sx["dec"] < 0]

icrf3k_nor = icrf3k[icrf3k["dec"] >= 0]
icrf3k_sou = icrf3k[icrf3k["dec"] < 0]

icrf3xka_nor = icrf3xka[icrf3xka["dec"] >= 0]
icrf3xka_sou = icrf3xka[icrf3xka["dec"] < 0]

gdr2_nor = gdr2[gdr2["dec"] >= 0]
gdr2_sou = gdr2[gdr2["dec"] < 0]

fig, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(
    figsize=(10, 12), ncols=2, nrows=3)

# Gaia DR2
ax0.plot(gdr2_sou["astrometric_n_obs_al"], gdr2_sou["pos_err"], "b*", ms=1)
ax0.plot(gdr2_nor["astrometric_n_obs_al"], gdr2_nor["pos_err"], "rx", ms=1)

ax0.set_yscale("log")
ax0.set_xscale("log")
ax0.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ Gaia DR2 ($\\mathrm{\\mu as}$)", fontsize=12)
ax0.set_xlabel("astrometric_n_obs_al", fontsize=12)

# Gaussian-like formal
y = 1 / np.sqrt(gdr2["astrometric_n_obs_al"]) * 1e3
ax0.plot(gdr2["astrometric_n_obs_al"], y, "k-", label="", lw=0.1)

ax0.set_ylim([10, 10000])
ax0.yaxis.set_ticks_position("both")
ax0.xaxis.set_ticks_position("both")

# ICRF1
ax1.plot(icrf1_nor["nb_del"], icrf1_nor["pos_err"], "rx", label="North", ms=1)
ax1.plot(icrf1_sou["nb_del"], icrf1_sou["pos_err"], "b*", label="South", ms=1)

# Gaussian-like formal
y = 3 / np.sqrt(icrf1["nb_del"]) * 1e3
ax1.plot(icrf1["nb_del"], y, "k-", label="", lw=0.1)

ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_ylim([100, 2000000])
ax1.yaxis.set_ticks_position("both")
ax1.xaxis.set_ticks_position("both")
ax1.legend(fontsize="x-small")
# ax1.text(100, 800, "ICRF1")
ax1.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF1 ($\\mathrm{\\mu as}$)", fontsize=12)
ax1.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF2
ax2.plot(icrf2_nor["nb_del"], icrf2_nor["pos_err"], "rx", label="North", ms=1)
ax2.plot(icrf2_sou["nb_del"], icrf2_sou["pos_err"], "b*", label="South", ms=1)

# Gaussian-like formal
y = 2 / np.sqrt(icrf2["nb_del"]) * 1e3
ax2.plot(icrf2["nb_del"], y, "k-", label="", lw=0.1)

ax2.set_yscale("log")
ax2.set_xscale("log")
ax2.set_ylim([20, 500000])
ax2.yaxis.set_ticks_position("both")
ax2.xaxis.set_ticks_position("both")
ax2.legend(fontsize="x-small")

ax2.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF2 ($\\mathrm{\\mu as}$)", fontsize=12)
ax2.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF3 SX
ax3.plot(icrf3sx_nor["nb_del"], icrf3sx_nor["pos_err"], "rx", label="North", ms=1)
ax3.plot(icrf3sx_sou["nb_del"], icrf3sx_sou["pos_err"], "b*", label="South", ms=1)

# Gaussian-like formal
y = 2 / np.sqrt(icrf3sx["nb_del"]) * 1e3
ax3.plot(icrf3sx["nb_del"], y, "k-", label="", lw=0.1)

ax3.set_ylim([20, 500000])

ax3.set_yscale("log")
ax3.set_xscale("log")
ax3.yaxis.set_ticks_position("both")
ax3.xaxis.set_ticks_position("both")
ax3.legend(fontsize="x-small")

ax3.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF3 S/X ($\\mathrm{\\mu as}$)", fontsize=12)
ax3.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF3 K
ax4.plot(icrf3k_nor["nb_del"], icrf3k_nor["pos_err"], "rx", label="North", ms=1)
ax4.plot(icrf3k_sou["nb_del"], icrf3k_sou["pos_err"], "b*", label="South", ms=1)

# Gaussian-like formal
y = 2 / np.sqrt(icrf3k["nb_del"]) * 1e3
ax4.plot(icrf3k["nb_del"], y, "k-", label="", lw=0.1)

ax4.set_ylim([40, 30000])
ax4.set_yscale("log")
ax4.set_xscale("log")
ax4.yaxis.set_ticks_position("both")
ax4.xaxis.set_ticks_position("both")
ax4.legend(fontsize="x-small")

ax4.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF3 K ($\\mathrm{\\mu as}$)", fontsize=12)
ax4.set_xlabel("$N_{\\rm delay}$", fontsize=12)

# ICRF3 XKa
ax5.plot(icrf3xka_nor["nb_del"],
         icrf3xka_nor["pos_err"], "rx", label="North", ms=1)
ax5.plot(icrf3xka_sou["nb_del"], icrf3xka_sou["pos_err"], "b*", label="South", ms=1)


# Gaussian-like formal
y = 1 / np.sqrt(icrf3xka["nb_del"]) * 1e3
ax5.plot(icrf3xka["nb_del"], y, "k-", label="", lw=0.1)

ax5.set_ylim([30, 200000])
ax5.set_yscale("log")
ax5.set_xscale("log")
ax5.yaxis.set_ticks_position("both")
ax5.xaxis.set_ticks_position("both")
ax5.legend(fontsize="x-small")

ax5.set_ylabel(
    "$\\sigma_{\\rm pos,max}$ ICRF3 X/Ka ($\\mathrm{\\mu as}$)", fontsize=12)
ax5.set_xlabel("$N_{\\rm delay}$", fontsize=12)

plt.tight_layout()
plt.savefig("../plots/error-vs-nobs-sou-nor_bgt.eps", hbox="tight")
# plt.savefig("../plots/fig02.eps")
# --------------------------------- END --------------------------------
