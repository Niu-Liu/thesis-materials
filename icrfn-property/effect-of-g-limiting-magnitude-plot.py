#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: effect-of-g-limiting-magnitude-plot.py
"""
Created on Sat Nov 16 09:51:46 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Some plots of the VSH parameters as a function of G-magnitude
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
minorLocator = MultipleLocator(10)
majorLocator = MultipleLocator(20)

# -----------------------------  MAINS -----------------------------
# Number of sources as a function of G-magnitude
gmag1, N0_1, N1_1, X0_1, rate1 = np.genfromtxt(
    "../logs/g-limiting-mag-stat-icrf1.log", unpack=True)
gmag2, N0_2, N1_2, X0_2, rate2 = np.genfromtxt(
    "../logs/g-limiting-mag-stat-icrf2.log", unpack=True)
gmag3sx, N0_3sx, N1_3sx, X0_3sx, rate3sx = np.genfromtxt(
    "../logs/g-limiting-mag-stat-icrf3sx.log", unpack=True)
gmag3k, N0_3k, N1_3k, X0_3k, rate3k = np.genfromtxt(
    "../logs/g-limiting-mag-stat-icrf3k.log", unpack=True)
gmag3xka, N0_3xka, N1_3xka, X0_3xka, rate3xka = np.genfromtxt(
    "../logs/g-limiting-mag-stat-icrf3xka.log", unpack=True)

fig, (ax0, ax1) = plt.subplots(figsize=(6, 4), nrows=2, sharex=True, sharey=True)

# Number of all common sources
ax0.plot(gmag1, N0_1, "g-v", label="ICRF1", ms=2, lw=0.5)
ax0.plot(gmag2, N0_2, "m-s", label="ICRF2", ms=2, lw=0.5)
ax0.plot(gmag3sx, N0_3sx, "r-x", label="ICRF3 S/X", ms=2, lw=0.5)
ax0.plot(gmag3k, N0_3k, "y-o", label="ICRF3 K", ms=2, lw=0.5)
ax0.plot(gmag3xka, N0_3xka, "b-*", label="ICRF3 X/Ka", ms=2, lw=0.5)

# ax0.axis([17.4, 21.1, 0, 3000])
# ax0.set_xlabel("$G_{\\rm max}$ (mag)")
ax0.set_ylabel("NO. common sources")
ax0.set_yscale("log")
# ax0.axis([17.4, 21.1, 100, 3000])

ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")

# ax0.grid()
ax0.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")

# Number of used sources
ax1.plot(gmag1, N1_1, "g-v", label="ICRF1", ms=2, lw=0.5)
ax1.plot(gmag2, N1_2, "m-s", label="ICRF2", ms=2, lw=0.5)
ax1.plot(gmag3sx, N1_3sx, "r-x", label="ICRF3 S/X", ms=2, lw=0.5)
ax1.plot(gmag3k, N1_3k, "y-o", label="ICRF3 K", ms=2, lw=0.5)
ax1.plot(gmag3xka, N1_3xka, "b-*", label="ICRF3 X/Ka", ms=2, lw=0.5)

# ax1.axis([17.4, 21.1, 0, 2500])
ax1.set_xlabel("$G_{\\rm max}$ (mag)")
ax1.set_ylabel("NO. used sources")
ax1.set_yscale("log")
ax0.axis([17.4, 21.1, 40, 3500])
ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")
# ax1.grid()

plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
# plt.savefig("../plots/nb-sou-on-g.eps")
plt.savefig("../plots/fig09.eps")
plt.close()

# VSH parameters
par1 = np.genfromtxt("../logs/g-limiting-mag-vsh02-icrf1.log")
par2 = np.genfromtxt("../logs/g-limiting-mag-vsh02-icrf2.log")
par3sx = np.genfromtxt("../logs/g-limiting-mag-vsh02-icrf3sx.log")
par3k = np.genfromtxt("../logs/g-limiting-mag-vsh02-icrf3k.log")
par3xka = np.genfromtxt("../logs/g-limiting-mag-vsh02-icrf3xka.log")

# --------------------- ROTATION -------------------------------
fig, (ax0, ax1, ax2) = plt.subplots(figsize=(6, 8), nrows=3, sharex=True)

# R1
ax0.errorbar(par1[:, 0], par1[:, 4], yerr=par1[:, 20],
             fmt="g-v", label="ICRF1", capsize=2, elinewidth=0.5, ms=3, lw=0.5)
ax0.errorbar(par2[:, 0], par2[:, 4], yerr=par2[:, 20],
             fmt="m-s", label="ICRF2", capsize=2, elinewidth=0.5, ms=3, lw=0.5)
ax0.errorbar(par3sx[:, 0], par3sx[:, 4], yerr=par3sx[:, 20],
             fmt="r-x", label="ICRF3 S/X", capsize=2, elinewidth=0.5, ms=3, lw=0.5)
ax0.errorbar(par3k[:, 0], par3k[:, 4], yerr=par3k[:, 20],
             fmt="y-o", label="ICRF3 K", capsize=2, elinewidth=0.5, ms=3, lw=0.5)
ax0.errorbar(par3xka[:, 0], par3xka[:, 4],
             yerr=par3xka[:, 20], fmt="b-*", label="ICRF3 X/Ka", capsize=2, elinewidth=0.5, ms=3, lw=0.5)

ax0.axis([17.4, 21.1, -80, 100])
ax0.set_ylabel("$R_1$ ($\mathrm{\mu}$as)")
ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")
ax0.yaxis.set_major_locator(majorLocator)
ax0.yaxis.set_minor_locator(minorLocator)
# ax0.grid()

# R2
ax1.errorbar(par1[:, 0], par1[:, 5], yerr=par1[:, 21],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par2[:, 0], par2[:, 5], yerr=par2[:, 21],
             fmt="m-s", label="ICRF2", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3sx[:, 0], par3sx[:, 5], yerr=par3sx[:, 21],
             fmt="r-x", label="ICRF3 S/X", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3k[:, 0], par3k[:, 5], yerr=par3k[:, 21],
             fmt="y-o", label="ICRF3 K", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3xka[:, 0], par3xka[:, 5],
             yerr=par3xka[:, 21], fmt="b-*", label="ICRF3 X/Ka", capsize=2, lw=0.5, elinewidth=0.5, ms=3)

ax1.axis([17.4, 21.1, -75, 50])
ax1.set_ylabel("$R_2$ ($\mathrm{\mu}$as)")
ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")
ax1.yaxis.set_minor_locator(minorLocator)
# ax1.grid()

# R3
ax2.errorbar(par1[:, 0], par1[:, 6], yerr=par1[:, 22],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par2[:, 0], par2[:, 6], yerr=par2[:, 22],
             fmt="m-s", label="ICRF2", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par3sx[:, 0], par3sx[:, 6], yerr=par3sx[:, 22],
             fmt="r-x", label="ICRF3 S/X", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par3k[:, 0], par3k[:, 6], yerr=par3k[:, 22],
             fmt="y-o", label="ICRF3 K", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par3xka[:, 0], par3xka[:, 6],
             yerr=par3xka[:, 22], fmt="b-*", label="ICRF3 X/Ka", capsize=2, lw=0.5, elinewidth=0.5, ms=3)

ax2.axis([17.4, 21.1, -75, 50])
ax2.set_xlabel("$G_{\\rm max}$ (mag)")
ax2.set_ylabel("$R_3$ ($\mathrm{\mu}$as)")
ax2.xaxis.set_ticks_position("both")
ax2.yaxis.set_ticks_position("both")
ax2.yaxis.set_minor_locator(minorLocator)
# ax2.grid()

ax0.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")
ax1.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")
ax2.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")

plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
# plt.savefig("../plots/rotation-on-g.eps")
plt.savefig("../plots/fig10.eps")
plt.close()

# --------------------- GLIDE -------------------------------
fig, (ax0, ax1, ax2) = plt.subplots(figsize=(6, 8), nrows=3, sharex=True)

# D1
ax0.errorbar(par1[:, 0], par1[:, 1], yerr=par1[:, 17],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par2[:, 0], par2[:, 1], yerr=par2[:, 17],
             fmt="m-s", label="ICRF2", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par3sx[:, 0], par3sx[:, 1], yerr=par3sx[:, 17],
             fmt="r-x", label="ICRF3 S/X", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par3k[:, 0], par3k[:, 1], yerr=par3k[:, 17],
             fmt="y-o", label="ICRF3 K", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par3xka[:, 0], par3xka[:, 1],
             yerr=par3xka[:, 20], fmt="b-*", label="ICRF3 X/Ka", capsize=2, lw=0.5, elinewidth=0.5, ms=3)

ax0.axis([17.4, 21.1, -110, 30])
ax0.set_ylabel("$D_1$ ($\mathrm{\mu}$as)")
ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")
# ax0.grid()

# D2
ax1.errorbar(par1[:, 0], par1[:, 2], yerr=par1[:, 18],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par2[:, 0], par2[:, 2], yerr=par2[:, 18],
             fmt="m-s", label="ICRF2", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3sx[:, 0], par3sx[:, 2], yerr=par3sx[:, 18],
             fmt="r-x", label="ICRF3 S/X", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3k[:, 0], par3k[:, 2], yerr=par3k[:, 18],
             fmt="y-o", label="ICRF3 K", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3xka[:, 0], par3xka[:, 2],
             yerr=par3xka[:, 18], fmt="b-*", label="ICRF3 X/Ka", capsize=2, lw=0.5, elinewidth=0.5, ms=3)

ax1.axis([17.4, 21.1, -40, 200])
ax1.set_ylabel("$D_2$ ($\mathrm{\mu}$as)")
ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")
# ax1.grid()

# D3
ax2.errorbar(par1[:, 0], par1[:, 3], yerr=par1[:, 19],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par2[:, 0], par2[:, 3], yerr=par2[:, 19],
             fmt="m-s", label="ICRF2", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par3sx[:, 0], par3sx[:, 3], yerr=par3sx[:, 19],
             fmt="r-x", label="ICRF3 S/X", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par3k[:, 0], par3k[:, 3], yerr=par3k[:, 19],
             fmt="y-o", label="ICRF3 K", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax2.errorbar(par3xka[:, 0], par3xka[:, 3],
             yerr=par3xka[:, 19], fmt="b-*", label="ICRF3 X/Ka", capsize=2, lw=0.5, elinewidth=0.5, ms=3)

ax2.axis([17.4, 21.1, -250, 330])
ax2.set_xlabel("$G_{\\rm max}$ (mag)")
ax2.set_ylabel("$D_3$ ($\mathrm{\mu}$as)")

ax2.xaxis.set_ticks_position("both")
ax2.yaxis.set_ticks_position("both")
# ax2.grid()

ax0.yaxis.set_minor_locator(minorLocator)
ax1.yaxis.set_minor_locator(minorLocator)
ax2.yaxis.set_minor_locator(MultipleLocator(20))

ax0.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")
ax1.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")
ax2.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")

plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
# plt.savefig("../plots/glide-on-g.eps")
plt.savefig("../plots/fig11.eps")
plt.close()

# E20 and M20
# E20
fig, (ax0, ax1) = plt.subplots(figsize=(6, 6), nrows=2, sharex=True)

ax0.errorbar(par1[:, 0], par1[:, 11], yerr=par1[:, 27],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par2[:, 0], par2[:, 11], yerr=par2[:, 27],
             fmt="m-s", label="ICRF2", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par3sx[:, 0], par3sx[:, 11], yerr=par3sx[:, 27],
             fmt="r-x", label="ICRF3 S/X", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par3k[:, 0], par3k[:, 11], yerr=par3k[:, 27],
             fmt="y-o", label="ICRF3 K", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par3xka[:, 0], par3xka[:, 11],
             yerr=par3xka[:, 27], fmt="b-*", label="ICRF3 X/Ka", capsize=2, lw=0.5, elinewidth=0.5, ms=3)

ax0.set_xlim([17.4, 21.1])
ax0.set_ylim([-140, 210])
# ax0.set_xlabel("$G_{\\rm max}$ (mag)")
ax0.set_ylabel("$E_{20}$ ($\mathrm{\mu}$as)")

ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")
ax0.yaxis.set_minor_locator(minorLocator)

# ax0.grid()
ax0.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")

# M20
ax1.errorbar(par1[:, 0], par1[:, 16], yerr=par1[:, 32],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par2[:, 0], par2[:, 16], yerr=par2[:, 32],
             fmt="m-s", label="ICRF2", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3sx[:, 0], par3sx[:, 16], yerr=par3sx[:, 32],
             fmt="r-x", label="ICRF3 S/X", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3k[:, 0], par3k[:, 16], yerr=par3k[:, 32],
             fmt="y-o", label="ICRF3 K", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax1.errorbar(par3xka[:, 0], par3xka[:, 16],
             yerr=par3xka[:, 32], fmt="b-*", label="ICRF3 X/Ka", capsize=2, lw=0.5, elinewidth=0.5, ms=3)

ax1.axis([17.4, 21.1, -70, 200])
ax1.set_xlabel("$G_{\\rm max}$ (mag)")
ax1.set_ylabel("$M_{20}$ ($\mathrm{\mu}$as)")

ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")

ax1.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")
ax1.yaxis.set_minor_locator(minorLocator)
# ax1.grid()

plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
# plt.savefig("../plots/EM20-on-g.eps")
plt.savefig("../plots/fig12.eps")
plt.close()

# E^{I,R}_2{1,2}
fig, (ax0, ax1, ax2, ax3) = plt.subplots(figsize=(6, 12), nrows=4, sharex=True)

# ER22
ax0.errorbar(par1[:, 0], par1[:, 7], yerr=par1[:, 23],
             fmt="g-v", label="ICRF1", capsize=2, lw=0.5, elinewidth=0.5, ms=3)
ax0.errorbar(par2[:, 0], par2[:, 7], yerr=par2[:, 23],
             fmt="m-s", label="ICRF2", capsize=2)
ax0.errorbar(par3sx[:, 0], par3sx[:, 7], yerr=par3sx[:, 23],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax0.errorbar(par3k[:, 0], par3k[:, 7], yerr=par3k[:, 23],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax0.errorbar(par3xka[:, 0], par3xka[:, 7],
             yerr=par3xka[:, 23], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax0.axis([17.4, 21.1, -50, 30])
ax0.set_xlabel("$G_{\\rm max}$ (mag)")
ax0.set_ylabel("$E^R_{22}$ ($\mathrm{\mu}$as)")

ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")

# ax0.grid()
# ax0.legend(fontsize="x-small")

# EI22
ax1.errorbar(par1[:, 0], par1[:, 8], yerr=par1[:, 24],
             fmt="g-v", label="ICRF1", capsize=2)
ax1.errorbar(par2[:, 0], par2[:, 8], yerr=par2[:, 24],
             fmt="m-s", label="ICRF2", capsize=2)
ax1.errorbar(par3sx[:, 0], par3sx[:, 8], yerr=par3sx[:, 24],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax1.errorbar(par3k[:, 0], par3k[:, 8], yerr=par3k[:, 24],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax1.errorbar(par3xka[:, 0], par3xka[:, 8],
             yerr=par3xka[:, 24], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax1.axis([17.4, 21.1, -30, 50])
ax1.set_xlabel("$G_{\\rm max}$ (mag)")
ax1.set_ylabel("$E^I_{22}$ ($\mathrm{\mu}$as)")

ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")

# ax1.grid()
# ax1.legend(fontsize="x-small")

# ER21
ax2.errorbar(par1[:, 0], par1[:, 9], yerr=par1[:, 25],
             fmt="g-v", label="ICRF1", capsize=2)
ax2.errorbar(par2[:, 0], par2[:, 9], yerr=par2[:, 25],
             fmt="m-s", label="ICRF2", capsize=2)
ax2.errorbar(par3sx[:, 0], par3sx[:, 9], yerr=par3sx[:, 25],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax2.errorbar(par3k[:, 0], par3k[:, 9], yerr=par3k[:, 25],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax2.errorbar(par3xka[:, 0], par3xka[:, 9],
             yerr=par3xka[:, 25], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax2.axis([17.4, 21.1, -150, 50])
ax2.set_xlabel("$G_{\\rm max}$ (mag)")
ax2.set_ylabel("$E^R_{21}$ ($\mathrm{\mu}$as)")

ax2.xaxis.set_ticks_position("both")
ax2.yaxis.set_ticks_position("both")

ax2.grid()
ax2.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")

# EI21
ax3.errorbar(par1[:, 0], par1[:, 10], yerr=par1[:, 26],
             fmt="g-v", label="ICRF1", capsize=2)
ax3.errorbar(par2[:, 0], par2[:, 10], yerr=par2[:, 26],
             fmt="m-s", label="ICRF2", capsize=2)
ax3.errorbar(par3sx[:, 0], par3sx[:, 10], yerr=par3sx[:, 26],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax3.errorbar(par3k[:, 0], par3k[:, 10], yerr=par3k[:, 26],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax3.errorbar(par3xka[:, 0], par3xka[:, 10],
             yerr=par3xka[:, 26], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax3.axis([17.4, 21.1, -120, 120])
ax3.set_xlabel("$G_{\\rm max}$ (mag)")
ax3.set_ylabel("$E^I_{21}$ ($\mathrm{\mu}$as)")

ax3.xaxis.set_ticks_position("both")
ax3.yaxis.set_ticks_position("both")

# ax3.grid()
# ax3.legend(fontsize="x-small")

plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
plt.savefig("../plots/vsh02-E-on-g.eps")
plt.close()

# MM^{I,R}_2{1}
fig, (ax0, ax1, ax2, ax3) = plt.subplots(figsize=(6, 12), nrows=4, sharex=True)

# MR22
ax0.errorbar(par1[:, 0], par1[:, 12], yerr=par1[:, 28],
             fmt="g-v", label="ICRF1", capsize=2)
ax0.errorbar(par2[:, 0], par2[:, 12], yerr=par2[:, 28],
             fmt="m-s", label="ICRF2", capsize=2)
ax0.errorbar(par3sx[:, 0], par3sx[:, 12], yerr=par3sx[:, 28],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax0.errorbar(par3k[:, 0], par3k[:, 12], yerr=par3k[:, 28],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax0.errorbar(par3xka[:, 0], par3xka[:, 12],
             yerr=par3xka[:, 28], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax0.axis([17.4, 21.1, -30, 70])
ax0.set_xlabel("$G_{\\rm max}$ (mag)")
ax0.set_ylabel("$M^R_{22}$ ($\mathrm{\mu}$as)")

ax0.xaxis.set_ticks_position("both")
ax0.yaxis.set_ticks_position("both")

ax0.grid()
# ax0.legend(fontsize="x-small")

# MI22
ax1.errorbar(par1[:, 0], par1[:, 13], yerr=par1[:, 29],
             fmt="g-v", label="ICRF1", capsize=2)
ax1.errorbar(par2[:, 0], par2[:, 13], yerr=par2[:, 29],
             fmt="m-s", label="ICRF2", capsize=2)
ax1.errorbar(par3sx[:, 0], par3sx[:, 13], yerr=par3sx[:, 29],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax1.errorbar(par3k[:, 0], par3k[:, 13], yerr=par3k[:, 29],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax1.errorbar(par3xka[:, 0], par3xka[:, 13],
             yerr=par3xka[:, 29], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax1.axis([17.4, 21.1, -50, 50])
ax1.set_xlabel("$G_{\\rm max}$ (mag)")
ax1.set_ylabel("$M^I_{22}$ ($\mathrm{\mu}$as)")

ax1.xaxis.set_ticks_position("both")
ax1.yaxis.set_ticks_position("both")

# ax1.grid()
# ax1.legend(fontsize="x-small")

# MR21
ax2.errorbar(par1[:, 0], par1[:, 14], yerr=par1[:, 30],
             fmt="g-v", label="ICRF1", capsize=2)
ax2.errorbar(par2[:, 0], par2[:, 14], yerr=par2[:, 30],
             fmt="m-s", label="ICRF2", capsize=2)
ax2.errorbar(par3sx[:, 0], par3sx[:, 14], yerr=par3sx[:, 30],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax2.errorbar(par3k[:, 0], par3k[:, 14], yerr=par3k[:, 30],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax2.errorbar(par3xka[:, 0], par3xka[:, 14],
             yerr=par3xka[:, 30], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax2.axis([17.4, 21.1, -50, 120])
ax2.set_xlabel("$G_{\\rm max}$ (mag)")
ax2.set_ylabel("$M^R_{21}$ ($\mathrm{\mu}$as)")

ax2.xaxis.set_ticks_position("both")
ax2.yaxis.set_ticks_position("both")

# ax2.grid()
ax2.legend(loc=2, bbox_to_anchor=(1.01, 1.0), borderaxespad=0., fontsize="xx-small")

# MI21
ax3.errorbar(par1[:, 0], par1[:, 15], yerr=par1[:, 31],
             fmt="g-v", label="ICRF1", capsize=2)
ax3.errorbar(par2[:, 0], par2[:, 15], yerr=par2[:, 31],
             fmt="m-s", label="ICRF2", capsize=2)
ax3.errorbar(par3sx[:, 0], par3sx[:, 15], yerr=par3sx[:, 31],
             fmt="r-x", label="ICRF3 S/X", capsize=2)
ax3.errorbar(par3k[:, 0], par3k[:, 15], yerr=par3k[:, 31],
             fmt="y-o", label="ICRF3 K", capsize=2)
ax3.errorbar(par3xka[:, 0], par3xka[:, 15],
             yerr=par3xka[:, 31], fmt="b-*", label="ICRF3 X/Ka", capsize=2)

ax3.axis([17.4, 21.1, -75, 100])
ax3.set_xlabel("$G_{\\rm max}$ (mag)")
ax3.set_ylabel("$M^I_{21}$ ($\mathrm{\mu}$as)")

ax3.xaxis.set_ticks_position("both")
ax3.yaxis.set_ticks_position("both")

# ax3.grid()
# ax3.legend(fontsize="x-small")

plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
plt.savefig("../plots/vsh02-M-on-g.eps")
plt.close()
# --------------------------------- END --------------------------------
