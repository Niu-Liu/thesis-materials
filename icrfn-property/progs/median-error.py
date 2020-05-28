#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: median-error.py
"""
Created on Sat May 18 11:07:37 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from matplotlib import pyplot as plt
import numpy as np

# My modules
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3
from my_progs.catalog.read_gaia import read_dr2_iers


# -----------------------------  FUNCTIONS -----------------------------
def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.00*height,
                "%d" % int(height), fontsize=8,
                ha="center", va="bottom")


# -----------------------------  MAINS -----------------------------
# ICRF1 catalog
icrf1 = read_icrf1()

# Fill the masked column
ra_err = icrf1["ra_err"]
dec_err = icrf1["dec_err"]
pos_err = icrf1["pos_err"]

# Median error
me_ra_i1 = np.median(ra_err[ra_err != 0]) * 1.e3
me_dec_i1 = np.median(dec_err[dec_err != 0]) * 1.e3
me_pos_i1 = np.median(pos_err[pos_err != 0]) * 1.e3

# ICRF2 catalog
icrf2 = read_icrf2()

# Median error
me_ra_i2 = np.median(icrf2["ra_err"]) * 1.e3
me_dec_i2 = np.median(icrf2["dec_err"]) * 1.e3
me_pos_i2 = np.median(icrf2["pos_err"]) * 1.e3

# ICRF3 S/X catalog
icrf3sx = read_icrf3(wv="sx")

# Median error
me_ra_i3sx = np.median(icrf3sx["ra_err"]) * 1.e3
me_dec_i3sx = np.median(icrf3sx["dec_err"]) * 1.e3
me_pos_i3sx = np.median(icrf3sx["pos_err"]) * 1.e3

# ICRF3 K catalog
icrf3k = read_icrf3(wv="k")

# Median error
me_ra_i3k = np.median(icrf3k["ra_err"]) * 1.e3
me_dec_i3k = np.median(icrf3k["dec_err"]) * 1.e3
me_pos_i3k = np.median(icrf3k["pos_err"]) * 1.e3

# ICRF3 X/Ka catalog
icrf3xka = read_icrf3(wv="xka")

# Median error
me_ra_i3xka = np.median(icrf3xka["ra_err"]) * 1.e3
me_dec_i3xka = np.median(icrf3xka["dec_err"]) * 1.e3
me_pos_i3xka = np.median(icrf3xka["pos_err"]) * 1.e3

# Gaia DR2 aux_iers catalog
gaiadr2 = read_dr2_iers()
# Bright sample
gaiadr2 = gaiadr2[gaiadr2["phot_g_mean_mag"] < 18.7]

# Median error
me_ra_g2 = np.median(gaiadr2["ra_err"]) * 1.e3
me_dec_g2 = np.median(gaiadr2["dec_err"]) * 1.e3
me_pos_g2 = np.median(gaiadr2["pos_err"]) * 1.e3

# Plot for median error
fig, ax = plt.subplots()
barwidth = 0.25

# data
catalogs = ["GAIA DR2", "ICRF1", "ICRF2", "ICRF3 S/X", "ICRF3 K", "ICRF3 X/KA"]
ra_pos = np.arange(len(catalogs)) - barwidth
dec_pos = np.arange(len(catalogs))
max_pos = np.arange(len(catalogs)) + barwidth

me_ra = np.array([me_ra_g2, me_ra_i1, me_ra_i2,
                  me_ra_i3sx, me_ra_i3k, me_ra_i3xka])
me_dec = np.array([me_dec_g2, me_dec_i1, me_dec_i2,
                   me_dec_i3sx, me_dec_i3k, me_dec_i3xka])
me_pos = np.array([me_pos_g2, me_pos_i1, me_pos_i2,
                   me_pos_i3sx, me_pos_i3k, me_pos_i3xka])

ra = ax.bar(ra_pos, me_ra,  width=barwidth, align="center",
            color="g", ecolor="black", label="$\\sigma_{\\alpha^*}$")
dec = ax.bar(dec_pos, me_dec,  width=barwidth, align="center",
             color="b", ecolor="black", label="$\\sigma_{\\delta}$")
eema = ax.bar(max_pos, me_pos,  width=barwidth, align="center",
              color="r", ecolor="black", label="$\\sigma_{\\rm pos,max}$")

autolabel(ra)
autolabel(dec)
autolabel(eema)

ax.set_xticks(dec_pos)
ax.set_xticklabels(catalogs, rotation="vertical")
ax.set_ylim([0, 800])
ax.set_ylabel("Median error ($\\mathrm{\\mu as}$)", fontsize=15)
# ax.yaxis.grid()  # horizontal lines
ax.legend()

plt.tight_layout()
# plt.savefig("../plots/median-error.eps")
plt.savefig("../plots/median-error_bgt.eps", hbox="tight")
# plt.savefig("../plots/fig01.eps")
# --------------------------------- END --------------------------------
