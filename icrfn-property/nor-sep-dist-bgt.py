#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: nor-sep-dist.py

"""
Created on Tue May 14 20:24:11 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, Column
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy import stats

# My modules
from my_progs.catalog.pos_diff import radio_cat_diff_calc, nor_sep_calc
from my_progs.catalog.read_gaia import read_dr2_iers
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3
from my_progs.catalog.vsh_deg2_cor import vsh_deg02_fitting, residual_calc02


# -----------------------------  FUNCTIONS -----------------------------
def calc_apri_nor_sep(com_sou):
    """Calculate the normalized separation after adjustments of VSH parameters
    """
    # Transform columns into np.array
    dra = np.array(com_sou["dra"])
    ddec = np.array(com_sou["ddec"])
    dra_err = np.array(com_sou["dra_err"])
    ddec_err = np.array(com_sou["ddec_err"])
    dra_ddec_cov = np.array(com_sou["dra_ddec_cov"])
    ra_rad = np.array(com_sou["ra"].to(u.radian))
    dec_rad = np.array(com_sou["dec"].to(u.radian))

    # Now re-calculate the normalized difference
    sep = nor_sep_calc(dra, dra_err, ddec, ddec_err,
                       dra_ddec_cov/dra_err/ddec_err)
    X = np.array(sep[-1])

    return X


def calc_post_nor_sep(com_sou, w2):
    """Calculate the normalized separation after adjustments of VSH parameters
    """
    # Transform columns into np.array
    dra = np.array(com_sou["dra"])
    ddec = np.array(com_sou["ddec"])
    dra_err = np.array(com_sou["dra_err"])
    ddec_err = np.array(com_sou["ddec_err"])
    dra_ddec_cov = np.array(com_sou["dra_ddec_cov"])
    ra_rad = np.array(com_sou["ra"].to(u.radian))
    dec_rad = np.array(com_sou["dec"].to(u.radian))

    # Remove the systematics
    dra_ns, ddec_ns = residual_calc02(dra, ddec, ra_rad, dec_rad, w2/1.e3)

    # Now re-calculate the normalized difference
    sep = nor_sep_calc(dra_ns, dra_err, ddec_ns, ddec_err,
                       dra_ddec_cov/dra_err/ddec_err)
    X = np.array(sep[-1])

    return X


# -----------------------------  MAIN -----------------------------
# ======================================================================
# Read the catalog
# ======================================================================
# Read Gaia DR2 IERS quasar data
gdr2 = read_dr2_iers()
gdr2 = gdr2[gdr2["phot_g_mean_mag"] < 18.7].filled()

# Read ICRF1 catalog
icrf1 = read_icrf1()

# Read ICRF2 catalog
icrf2 = read_icrf2()

# Read ICRF3 S/X catalog
icrf3sx = read_icrf3(wv="sx")

# Read ICRF3 K catalog
icrf3k = read_icrf3(wv="k")

# Read ICRF3 S/X catalog
icrf3xka = read_icrf3(wv="xka")

# ======================================================================
# Calculate the position difference
# ======================================================================
# ICRF1 - Gaia DR2
icrf1_gdr2 = radio_cat_diff_calc(icrf1, gdr2, "iers_name", ["icrf1", "gdr2"])

# ICRF2 - Gaia DR2
icrf2_gdr2 = radio_cat_diff_calc(icrf2, gdr2, "iers_name", ["icrf2", "gdr2"])

# ICRF3SX - Gaia DR2
icrf3sx_gdr2 = radio_cat_diff_calc(
    icrf3sx, gdr2, "iers_name", ["icrf3sx", "gdr2"])

# ICRF3K - Gaia DR2
icrf3k_gdr2 = radio_cat_diff_calc(
    icrf3k, gdr2, "iers_name", ["icrf3k", "gdr2"])

# ICRF3XKa - Gaia DR2
icrf3xka_gdr2 = radio_cat_diff_calc(
    icrf3xka, gdr2, "iers_name", ["icrf3xka", "gdr2"])


# ======================================================================
# VSH parameters
# ======================================================================
w2_icrf1 = np.genfromtxt(
    "../logs/icrf1_gaiadr2_bgt_vsh02.log", usecols=(1,), skip_header=1)
w2_icrf2 = np.genfromtxt(
    "../logs/icrf2_gaiadr2_bgt_vsh02.log", usecols=(1,), skip_header=1)
w2_icrf3sx = np.genfromtxt(
    "../logs/icrf3sx_gaiadr2_bgt_vsh02.log", usecols=(1,), skip_header=1)
w2_icrf3k = np.genfromtxt(
    "../logs/icrf3k_gaiadr2_bgt_vsh02.log", usecols=(1,), skip_header=1)
w2_icrf3xka = np.genfromtxt(
    "../logs/icrf3xka_gaiadr2_bgt_vsh02.log", usecols=(1,), skip_header=1)

# ======================================================================
# Add adjustments
# ======================================================================
# Apriori fit
apri_x_icrf1 = calc_apri_nor_sep(icrf1_gdr2)
apri_x_icrf2 = calc_apri_nor_sep(icrf2_gdr2)
apri_x_icrf3sx = calc_apri_nor_sep(icrf3sx_gdr2)
apri_x_icrf3k = calc_apri_nor_sep(icrf3k_gdr2)
apri_x_icrf3xka = calc_apri_nor_sep(icrf3xka_gdr2)

# Post-fit
post_x_icrf1 = calc_post_nor_sep(icrf1_gdr2, w2_icrf1)
post_x_icrf2 = calc_post_nor_sep(icrf2_gdr2, w2_icrf2)
post_x_icrf3sx = calc_post_nor_sep(icrf3sx_gdr2, w2_icrf3sx)
post_x_icrf3k = calc_post_nor_sep(icrf3k_gdr2, w2_icrf3k)
post_x_icrf3xka = calc_post_nor_sep(icrf3xka_gdr2, w2_icrf3xka)

# ======================================================================
# Make the plot
# ======================================================================
# A (standard) Rayleigh distribution
bins_array = np.linspace(0, 10, 50)
rayleigh_dist = stats.rayleigh.pdf(bins_array) * 10. / 50 * 100

minorLocator = MultipleLocator(1)

# Plot the apriori/postori X
# Distribution of apriori-fit normalized separation
fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(figsize=(6, 10), nrows=5,
                                              sharex=True, sharey=True)

# ICRF1
weights = np.ones_like(apri_x_icrf1) * 100. / apri_x_icrf1.size
ax0.hist(apri_x_icrf1, bins_array, weights=weights,
         color="grey", histtype="stepfilled", label="pre-fit", alpha=0.5)
ax0.hist(post_x_icrf1, bins_array, weights=weights,
         color="r", histtype="step", label="post-fit")
ax0.plot(bins_array, rayleigh_dist, "b--", linewidth=1)
ax0.text(4, 12, "ICRF1", fontsize=15)
ax0.legend(fontsize="small")

ax0.set_ylim([0, 16])
ax0.set_xlim([0, 10])
ax0.set_ylabel("% in bin", fontsize=12)
# ax0.grid()
ax0.set_yticks(np.arange(0, 15, 5))
ax0.yaxis.set_minor_locator(minorLocator)

# ICRF2
weights = np.ones_like(apri_x_icrf2) * 100. / apri_x_icrf2.size
ax1.hist(apri_x_icrf2, bins_array, weights=weights,
         color="grey", histtype="stepfilled", label="pre-fit")
ax1.hist(post_x_icrf2, bins_array, weights=weights,
         color="r", histtype="step", label="post-fit")
ax1.plot(bins_array, rayleigh_dist, "b--", linewidth=1)
ax1.text(4, 12, "ICRF2", fontsize=15)
ax1.legend(fontsize="small")

# ax1.set_ylim([0, 13])
ax1.set_ylabel("% in bin", fontsize=12)
# ax1.grid()
ax1.set_yticks(np.arange(0, 15, 5))
ax1.yaxis.set_minor_locator(minorLocator)
ax1.xaxis.set_minor_locator(minorLocator)

# ICRF3 SX
weights = np.ones_like(apri_x_icrf3sx) * 100. / apri_x_icrf3sx.size
ax2.hist(apri_x_icrf3sx, bins_array, weights=weights,
         color="grey", histtype="stepfilled", label="pre-fit")
ax2.hist(post_x_icrf3sx, bins_array, weights=weights,
         color="r", histtype="step", label="post-fit")
ax2.plot(bins_array, rayleigh_dist, "b--", linewidth=1)
ax2.text(4, 12, "ICRF3 S/X", fontsize=15)
ax2.legend(fontsize="small")

# ax2.set_ylim([0, 13])
ax2.set_ylabel("% in bin", fontsize=12)
# ax2.grid()
ax2.set_yticks(np.arange(0, 15, 5))
ax2.yaxis.set_minor_locator(minorLocator)

# ICRF3 K
weights = np.ones_like(apri_x_icrf3k) * 100. / apri_x_icrf3k.size
ax3.hist(apri_x_icrf3k, bins_array, weights=weights,
         color="grey", histtype="stepfilled", label="pre-fit")
ax3.hist(post_x_icrf3k, bins_array, weights=weights,
         color="r", histtype="step", label="post-fit")
ax3.plot(bins_array, rayleigh_dist, "b--", linewidth=1)
ax3.text(4, 12, "ICRF3 K", fontsize=15)
ax3.legend(fontsize="small")

# ax3.set_ylim([0, 13])
ax3.set_ylabel("% in bin", fontsize=12)
# ax3.grid()
ax3.set_yticks(np.arange(0, 15, 5))
ax3.yaxis.set_minor_locator(minorLocator)

# ICRF3 X/Ka
weights = np.ones_like(apri_x_icrf3xka) * 100. / apri_x_icrf3xka.size
ax4.hist(apri_x_icrf3xka, bins_array, weights=weights,
         color="grey", histtype="stepfilled", label="pre-fit")
ax4.hist(post_x_icrf3xka, bins_array, weights=weights,
         color="r", histtype="step", label="post-fit")
ax4.plot(bins_array, rayleigh_dist, "b--", linewidth=1)
ax4.text(4, 12, "ICRF3 X/Ka", fontsize=15)
ax4.legend(fontsize="small")

# ax4.set_ylim([0, 13])
ax4.set_xlabel("Normalized seperation $X$", fontsize=12)
ax4.set_ylabel("% in bin", fontsize=12)
# ax4.grid()
ax4.set_yticks(np.arange(0, 15, 5))
ax4.yaxis.set_minor_locator(minorLocator)

plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig("../plots/nor-sep-dist_bgt.eps", hbox="tight")
plt.close()
# --------------------------------- END --------------------------------
