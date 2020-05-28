#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: vsh_param.py
"""
Created on Wed May 15 08:50:54 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from my_progs.catalog.vec_mod import vec_mod_calc


# -----------------------------  MAINS -----------------------------
# ======================================================================
# VSH parameters
# ======================================================================
# ICRF1
# param1, perr1 = np.genfromtxt(
# "../logs/icrf1_gaiadr2_vsh02.log",
# usecols=(1, 2), skip_header=1, unpack=True)
param1, perr1 = np.genfromtxt(
    "../logs/icrf1_gaiadr2_bgt_vsh02.log",
    usecols=(1, 2), skip_header=1, unpack=True)

# ICRF2
# param2, perr2 = np.genfromtxt(
#     "../logs/icrf2_gaiadr2_vsh02.log",
#     usecols=(1, 2), skip_header=1, unpack=True)
param2, perr2 = np.genfromtxt(
    "../logs/icrf2_gaiadr2_bgt_vsh02.log",
    usecols=(1, 2), skip_header=1, unpack=True)

# ICRF3 S/X
# param3sx, perr3sx = np.genfromtxt(
#     "../logs/icrf3sx_gaiadr2_vsh02.log",
#     usecols=(1, 2), skip_header=1, unpack=True)
param3sx, perr3sx = np.genfromtxt(
    "../logs/icrf3sx_gaiadr2_bgt_vsh02.log",
    usecols=(1, 2), skip_header=1, unpack=True)

# ICRF3 K
# param3k, perr3k = np.genfromtxt(
#     "../logs/icrf3k_gaiadr2_vsh02.log",
#     usecols=(1, 2), skip_header=1, unpack=True)
param3k, perr3k = np.genfromtxt(
    "../logs/icrf3k_gaiadr2_bgt_vsh02.log",
    usecols=(1, 2), skip_header=1, unpack=True)

# ICRF3 X/Ka
# param3xka, perr3xka = np.genfromtxt(
#     "../logs/icrf3xka_gaiadr2_vsh02.log",
#     usecols=(1, 2), skip_header=1, unpack=True)
param3xka, perr3xka = np.genfromtxt(
    "../logs/icrf3xka_gaiadr2_bgt_vsh02.log",
    usecols=(1, 2), skip_header=1, unpack=True)

# ======================================================================
# Separate the parameters into rotation, glide, and quadrupole terms
# ======================================================================
# ICRF1
gli1 = param1[:3]
rot1 = param1[3:6]
qua1 = param1[6:]

gerr1 = perr1[:3]
rerr1 = perr1[3:6]
qerr1 = perr1[6:]

glimod1, glierr1 = vec_mod_calc(gli1, gerr1)
rotmod1, roterr1 = vec_mod_calc(rot1, rerr1)

# ICRF2
gli2 = param2[:3]
rot2 = param2[3:6]
qua2 = param2[6:]

gerr2 = perr2[:3]
rerr2 = perr2[3:6]
qerr2 = perr2[6:]

glimod2, glierr2 = vec_mod_calc(gli2, gerr2)
rotmod2, roterr2 = vec_mod_calc(rot2, rerr2)


# ICRF3 SX
gli3sx = param3sx[:3]
rot3sx = param3sx[3:6]
qua3sx = param3sx[6:]

gerr3sx = perr3sx[:3]
rerr3sx = perr3sx[3:6]
qerr3sx = perr3sx[6:]

glimod3sx, glierr3sx = vec_mod_calc(gli3sx, gerr3sx)
rotmod3sx, roterr3sx = vec_mod_calc(rot3sx, rerr3sx)

# ICRF3 K
gli3k = param3k[:3]
rot3k = param3k[3:6]
qua3k = param3k[6:]

gerr3k = perr3k[:3]
rerr3k = perr3k[3:6]
qerr3k = perr3k[6:]

glimod3k, glierr3k = vec_mod_calc(gli3k, gerr3k)
rotmod3k, roterr3k = vec_mod_calc(rot3k, rerr3k)

# ICRF3 X/Ka
gli3xka = param3xka[:3]
rot3xka = param3xka[3:6]
qua3xka = param3xka[6:]

gerr3xka = perr3xka[:3]
rerr3xka = perr3xka[3:6]
qerr3xka = perr3xka[6:]

glimod3xka, glierr3xka = vec_mod_calc(gli3xka, gerr3xka)
rotmod3xka, roterr3xka = vec_mod_calc(rot3xka, rerr3xka)

# ======================================================================
# Plot
# ======================================================================
minorLocator = MultipleLocator(10)
majorLocator = MultipleLocator(50)

# Rotation
fig, ax = plt.subplots()

barwidth = 0.15
loc = 0.15

terms = ["$R_1$", "$R_2$", "$R_3$", "$R$"]

icrf1_pos = np.arange(len(terms)) - 2 * loc
icrf2_pos = np.arange(len(terms)) - 1 * loc
icrf3sx_pos = np.arange(len(terms))
icrf3k_pos = np.arange(len(terms)) + 1 * loc
icrf3xka_pos = np.arange(len(terms)) + 2 * loc

par1 = np.concatenate((rot1, [rotmod1]))
err1 = np.concatenate((rerr1, [roterr1]))

par2 = np.concatenate((rot2, [rotmod2]))
err2 = np.concatenate((rerr2, [roterr2]))

par3sx = np.concatenate((rot3sx, [rotmod3sx]))
err3sx = np.concatenate((rerr3sx, [roterr3sx]))

par3k = np.concatenate((rot3k, [rotmod3k]))
err3k = np.concatenate((rerr3k, [roterr3k]))

par3xka = np.concatenate((rot3xka, [rotmod3xka]))
err3xka = np.concatenate((rerr3xka, [roterr3xka]))

ax.bar(icrf1_pos, 2 * err1, bottom=par1-err1, width=barwidth,
       color="g", ecolor="black", label="ICRF1")
ax.bar(icrf2_pos, 2 * err2, bottom=par2-err2, width=barwidth,
       color="m", ecolor="black", label="ICRF2")
ax.bar(icrf3sx_pos, 2 * err3sx, bottom=par3sx-err3sx, width=barwidth,
       color="r", ecolor="black", label="ICRF3 S/X")
ax.bar(icrf3k_pos, 2 * err3k, bottom=par3k-err3k, width=barwidth,
       color="y", ecolor="black", label="ICRF3 K")
ax.bar(icrf3xka_pos, 2 * err3xka, bottom=par3xka-err3xka, width=barwidth,
       color="b", ecolor="black", label="ICRF3 X/Ka")

ax.hlines(par1, icrf1_pos-0.5*barwidth, icrf1_pos+0.5*barwidth, color="k")
ax.hlines(par2, icrf2_pos-0.5*barwidth, icrf2_pos+0.5*barwidth, color="k")
ax.hlines(par3sx, icrf3sx_pos-0.5*barwidth, icrf3sx_pos+0.5*barwidth, color="k")
ax.hlines(par3k, icrf3k_pos-0.5*barwidth, icrf3k_pos+0.5*barwidth, color="k")
ax.hlines(par3xka, icrf3xka_pos-0.5*barwidth, icrf3xka_pos+0.5*barwidth, color="k")

# ax.plot(icrf1_pos, par1, "kx")
# ax.plot(icrf2_pos, par2, "kx")
# ax.plot(icrf3sx_pos, par3sx, "kx")
# ax.plot(icrf3k_pos, par3k, "kx")
# ax.plot(icrf3xka_pos, par3xka, "kx")
ax.vlines(2.5, -70, 90, lw=1, linestyles="dashdot", colors="r")

ax.set_xticks(np.arange(len(terms)))
ax.set_xticklabels(terms, fontsize=15)
ax.set_xlim([-0.5, len(icrf1_pos)-0.5])
ax.set_ylim([-70, 90])
ax.set_ylabel("Rotation ($\\mathrm{\mu as}$)", fontsize=15)
ax.xaxis.set_ticks_position("both")
ax.yaxis.set_ticks_position("both")
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
# ax.xaxis.grid(True, which="Minor")
# ax.yaxis.grid()  # horizontal lines
ax.yaxis.set_minor_locator(minorLocator)
ax.legend(fontsize="x-small")

plt.tight_layout()
plt.savefig("../plots/rotation_bgt.eps", hbox="tight")
plt.close()

# Glide
fig, ax = plt.subplots()

# barwidth = 0.1
# loc = 0.12

terms = ["$D_1$", "$D_2$", "$D_3$", "$D$"]

icrf1_pos = np.arange(len(terms)) - 2 * loc
icrf2_pos = np.arange(len(terms)) - 1 * loc
icrf3sx_pos = np.arange(len(terms))
icrf3k_pos = np.arange(len(terms)) + 1 * loc
icrf3xka_pos = np.arange(len(terms)) + 2 * loc

par1 = np.concatenate((gli1, [glimod1]))
err1 = np.concatenate((gerr1, [glierr1]))

par2 = np.concatenate((gli2, [glimod2]))
err2 = np.concatenate((gerr2, [glierr2]))

par3sx = np.concatenate((gli3sx, [glimod3sx]))
err3sx = np.concatenate((gerr3sx, [glierr3sx]))

par3k = np.concatenate((gli3k, [glimod3k]))
err3k = np.concatenate((gerr3k, [glierr3k]))

par3xka = np.concatenate((gli3xka, [glimod3xka]))
err3xka = np.concatenate((gerr3xka, [glierr3xka]))

ax.bar(icrf1_pos, 2 * err1, bottom=par1-err1, width=barwidth,
       color="g", ecolor="black", label="ICRF1")
ax.bar(icrf2_pos, 2 * err2, bottom=par2-err2, width=barwidth,
       color="m", ecolor="black", label="ICRF2")
ax.bar(icrf3sx_pos, 2 * err3sx, bottom=par3sx-err3sx, width=barwidth,
       color="r", ecolor="black", label="ICRF3 S/X")
ax.bar(icrf3k_pos, 2 * err3k, bottom=par3k-err3k, width=barwidth,
       color="y", ecolor="black", label="ICRF3 K")
ax.bar(icrf3xka_pos, 2 * err3xka, bottom=par3xka-err3xka, width=barwidth,
       color="b", ecolor="black", label="ICRF3 X/Ka")

ax.hlines(par1, icrf1_pos-0.5*barwidth, icrf1_pos+0.5*barwidth, color="k")
ax.hlines(par2, icrf2_pos-0.5*barwidth, icrf2_pos+0.5*barwidth, color="k")
ax.hlines(par3sx, icrf3sx_pos-0.5*barwidth, icrf3sx_pos+0.5*barwidth, color="k")
ax.hlines(par3k, icrf3k_pos-0.5*barwidth, icrf3k_pos+0.5*barwidth, color="k")
ax.hlines(par3xka, icrf3xka_pos-0.5*barwidth, icrf3xka_pos+0.5*barwidth, color="k")

# ax.plot(icrf1_pos, par1, "kx")
# ax.plot(icrf2_pos, par2, "kx")
# ax.plot(icrf3sx_pos, par3sx, "kx")
# ax.plot(icrf3k_pos, par3k, "kx")
# ax.plot(icrf3xka_pos, par3xka, "kx")
ax.vlines(2.5, -250, 300, lw=1, linestyles="dashdot", colors="r")

ax.set_xticks(np.arange(len(terms)))
ax.set_xticklabels(terms, fontsize=15)
ax.set_xlim([-0.5, len(icrf1_pos)-0.5])
ax.set_ylim([-250, 300])
ax.set_ylabel("Glide ($\\mathrm{\mu as}$)", fontsize=15)
ax.xaxis.set_ticks_position("both")
ax.yaxis.set_ticks_position("both")
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
# ax.xaxis.grid(True, which="Minor")
# ax.yaxis.grid()  # horizontal lines
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator)
ax.legend(fontsize="x-small")

plt.tight_layout()
# plt.savefig("../plots/glide.eps")
plt.savefig("../plots/glide_bgt.eps", hbox="tight")
# plt.savefig("../plots/fig07.eps")
plt.close()

# Quadruple terms

fig, ax = plt.subplots()

# barwidth = 0.2
# loc = 0.22

terms = ["$E_{22}^R$", "$E_{22}^I$", "$E_{21}^R$", "$E_{21}^I$", "$E_{20}$",
         "$M_{22}^R$", "$M_{22}^I$", "$M_{21}^R$", "$M_{21}^I$", "$M_{20}$"]

icrf1_pos = np.arange(len(terms)) - 2 * loc
icrf2_pos = np.arange(len(terms)) - 1 * loc
icrf3sx_pos = np.arange(len(terms))
icrf3k_pos = np.arange(len(terms)) + 1 * loc
icrf3xka_pos = np.arange(len(terms)) + 2 * loc

ax.bar(icrf1_pos, 2 * qerr1, bottom=qua1-qerr1, width=barwidth,
       color="g", ecolor="black", label="ICRF1")
ax.bar(icrf2_pos, 2 * qerr2, bottom=qua2-qerr2, width=barwidth,
       color="m", ecolor="black", label="ICRF2")
ax.bar(icrf3sx_pos, 2 * qerr3sx, bottom=qua3sx-qerr3sx, width=barwidth,
       color="r", ecolor="black", label="ICRF3 S/X")
ax.bar(icrf3k_pos, 2 * qerr3k, bottom=qua3k-qerr3k, width=barwidth,
       color="y", ecolor="black", label="ICRF3 K")
ax.bar(icrf3xka_pos, 2 * qerr3xka, bottom=qua3xka-qerr3xka, width=barwidth,
       color="b", ecolor="black", label="ICRF3 X/Ka")

ax.hlines(qua1, icrf1_pos-0.5*barwidth, icrf1_pos+0.5*barwidth, color="k")
ax.hlines(qua2, icrf2_pos-0.5*barwidth, icrf2_pos+0.5*barwidth, color="k")
ax.hlines(qua3sx, icrf3sx_pos-0.5*barwidth, icrf3sx_pos+0.5*barwidth, color="k")
ax.hlines(qua3k, icrf3k_pos-0.5*barwidth, icrf3k_pos+0.5*barwidth, color="k")
ax.hlines(qua3xka, icrf3xka_pos-0.5*barwidth, icrf3xka_pos+0.5*barwidth, color="k")

# ax.plot(icrf1_pos, qua1, "kx")
# ax.plot(icrf2_pos, qua2, "kx")
# ax.plot(icrf3sx_pos, qua3sx, "kx")
# ax.plot(icrf3k_pos, qua3k, "kx")
# ax.plot(icrf3xka_pos, qua3xka, "kx")

ax.set_xticks(np.arange(len(terms)))
ax.set_xticklabels(terms, fontsize=15)
ax.set_xlim([-0.5, len(icrf1_pos)-0.5])
# ax.set_ylim([-100, 200])
ax.set_ylim([-120, 180])
ax.set_ylabel("Quadrupole ($\\mathrm{\mu as}$)", fontsize=15)
ax.xaxis.set_ticks_position("both")
ax.yaxis.set_ticks_position("both")
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
# ax.xaxis.grid(True, which="Minor")
# ax.yaxis.grid()  # horizontal lines
ax.yaxis.set_minor_locator(minorLocator)
ax.legend(fontsize="x-small")

plt.tight_layout()
# plt.savefig("../plots/quadrupole.eps")
plt.savefig("../plots/quadrupole_bgt.eps", hbox="tight")
# plt.savefig("../plots/fig08.eps")
plt.close()
# --------------------------------- END --------------------------------
