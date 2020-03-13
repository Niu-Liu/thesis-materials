# File name: angular-separation-icrf-vs-icrf3.py
"""
Created on Sat May 18 11:37:11 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""


from astropy.table import join, Table
import astropy.units as u
from matplotlib import pyplot as plt
import numpy as np

# My modules
from my_progs.catalog.read_icrf import read_icrf1, read_icrf2, read_icrf3
from my_progs.catalog.read_gaia import read_dr2_iers
from my_progs.catalog.pos_diff import radio_cat_diff_calc, nor_sep_calc
from my_progs.catalog.vsh_deg2_cor import residual_calc02


# -----------------------------  FUNCTIONS -----------------------------
# Define a function to calculate the angular/normalized separation.
def calc_apri_sep(com_sou):
    """Calculate the normalized separation after adjustments of VSH parameters
    """
    # Transform columns into np.array
    # Transform columns into np.array
    dra = np.array(com_sou["dra"])
    ddec = np.array(com_sou["ddec"])
    dra_err = np.array(com_sou["dra_err"])
    ddec_err = np.array(com_sou["ddec_err"])
    dra_ddec_cov = np.array(com_sou["dra_ddec_cov"])
    ra_rad = np.array(com_sou["ra"].to(u.radian))
    dec_rad = np.array(com_sou["dec"].to(u.radian))

    # Now re-calculate the normalized difference
    sep = nor_sep_calc(
        dra, dra_err, ddec, ddec_err, dra_ddec_cov/dra_err/ddec_err)

    sep_table = Table([com_sou["iers_name"], sep[0], sep[-1]],
                      names=["iers_name", "ang_sep", "nor_sep"])
    return sep_table


def calc_post_sep(com_sou, vshpar):
    """Calculate the angular/normalized separation after adjustments of VSH parameters
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
    dra_ns, ddec_ns = residual_calc02(dra, ddec, ra_rad, dec_rad, vshpar/1.e3)

    # Now re-calculate the normalized difference
    sep = nor_sep_calc(
        dra_ns, dra_err, ddec_ns, ddec_err, dra_ddec_cov/dra_err/ddec_err)

    sep_table = Table([com_sou["iers_name"], sep[0], sep[-1]],
                      names=["iers_name", "ang_sep", "nor_sep"])

    return sep_table


# ------------------------------   MAIN -- -----------------------------
# Read the catalogs.

# ICRF2 catalog
icrf2 = read_icrf2()

# ICRF3 S/X catalog
icrf3sx = read_icrf3(wv="sx")

# Gaia DR2 aux_iers catalog
gdr2 = read_dr2_iers()


# Calculate the positional differences.
# ICRF2 - Gaia DR2
icrf2_gdr2 = radio_cat_diff_calc(icrf2, gdr2, "iers_name")

# ICRF3 SX - Gaia DR2
icrf3sx_gdr2 = radio_cat_diff_calc(icrf3sx, gdr2, "iers_name")


# VSH parameters
vshpar_icrf2 = np.genfromtxt("../logs/icrf2_gaiadr2_vsh02.log",
                             usecols=(1, ),
                             skip_header=1)
vshpar_icrf3sx = np.genfromtxt("../logs/icrf3sx_gaiadr2_vsh02.log",
                               usecols=(1, ),
                               skip_header=1)


# Calculate the separation after removing the systematics.
post_sep_icrf2 = calc_post_sep(icrf2_gdr2, vshpar_icrf2)
post_sep_icrf3sx = calc_post_sep(icrf3sx_gdr2, vshpar_icrf3sx)

# Crossmatch the two tables.
com_sou_2v3sx = join(post_sep_icrf2, post_sep_icrf3sx,
                     keys="iers_name", table_names=["2", "3sx"])

# Unit transformation
com_sou_2v3sx["ang_sep_2"] = com_sou_2v3sx["ang_sep_2"] * 1.e3
com_sou_2v3sx["ang_sep_3sx"] = com_sou_2v3sx["ang_sep_3sx"] * 1.e3

# Plot the angular separation of ICRF3 vs. ICRF2.
fig, (ax0, ax1) = plt.subplots(figsize=(8, 4), ncols=2)

# Angular separation
ax0.plot(com_sou_2v3sx["ang_sep_2"], com_sou_2v3sx["ang_sep_3sx"], "b+", ms=1)
# Plot y=x
x = np.arange(10, 500000, 2**0.2)
ax0.plot(x, x, "k", lw=0.2)

ax0.set_xlim([10, 500000])
ax0.set_ylim([10, 500000])
ax0.set_xscale("log")
ax0.set_yscale("log")

ax0.set_xlabel("$\\rho_{\\rm ICRF2}$ ($\\mathrm{\\mu as}$)")
ax0.set_ylabel("$\\rho_{\\rm ICRF3\ S/X}$ ($\\mathrm{\\mu as}$)")

ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")

# ax0.grid()

# Normalized separation
ax1.plot(com_sou_2v3sx["nor_sep_2"], com_sou_2v3sx["nor_sep_3sx"], "b+", ms=1)
# Plot y=x
x = np.arange(0.01, 400, 2**0.2)
ax1.plot(x, x, "k", lw=0.2)

ax1.set_xlim([0.01, 400])
ax1.set_ylim([0.01, 400])
ax1.set_xscale("log")
ax1.set_yscale("log")

ax1.set_xlabel("$X_{\\rm ICRF2}$")
ax1.set_ylabel("$X_{\\rm ICRF3\ S/X}$")

ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")

# ax1.grid()

plt.tight_layout()
plt.savefig("../plots/post_sep_icrf2_vs_icrf3.eps", hbox="tight")
plt.close()
# --------------------------------- END --------------------------------
