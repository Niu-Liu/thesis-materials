#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: glide_plot.py
"""
Created on Wed Jun 13 10:40:04 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos


# -----------------------------  FUNCTIONS -----------------------------
def glide_field_gen(gv, ra, dec):
    """Generate a field at (ra,dec) the glide scale.

    Parameters
    ----------
    gv : array of float
        glide vector
    ra : array of float
        right ascension
    dec : array of float
        declination

    Returns
    ----------
    g_dra : array of float
        RA offset induced by glide
    g_ddec : array of float
        Dec. offset induced by glide
    """

    g1, g2, g3 = gv
    g_dra = -g1 * sin(ra) + g2 * cos(ra)
    g_ddec = - g1*cos(ra)*sin(dec) - g2*sin(ra)*sin(dec) + g3 * cos(dec)

    return g_dra, g_ddec


def glide_plot(gv, output=None, fig_title=None):
    """Plot fot glide field.

    Parameters
    ----------
    gv : array of float
        glide vector
    output : string
        Full name of the output file
    fig_title : string
        Title of the figure

    Returns
    ----------
    None
    """

    ra = np.linspace(0, 360, 20) * u.deg
    dec = (np.linspace(-90, 90, 20)) * u.deg
    c = SkyCoord(ra=ra, dec=dec, frame='icrs')

    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian

    # Glide field
    X, Y = np.meshgrid(ra_rad, dec_rad)
    U, V = glide_field_gen(gv, X, Y)

    # GC position
    c1_gal = SkyCoord(l=0 * u.deg, b=0 * u.deg, frame="galactic")
    c1_icrs = c1_gal.icrs
    ra_gc_rad = c1_icrs.ra.wrap_at(180 * u.deg).radian
    dec_gc_rad = c1_icrs.dec.radian

    # Anti-GC position
    c2_gal = SkyCoord(l=180 * u.deg, b=0 * u.deg, frame="galactic")
    c2_icrs = c2_gal.icrs
    ra_agc_rad = c2_icrs.ra.wrap_at(180 * u.deg).radian
    dec_agc_rad = c2_icrs.dec.radian

    # Plot
    plt.figure(figsize=(8, 4.2))
    plt.subplot(111, projection="aitoff")

    Q = plt.quiver(X, Y, U, V, units='xy', scale=100.0)
    qk = plt.quiverkey(Q, 0.90, 0.90, 50, r'$50 \mu as$', labelpos='E',
                       coordinates='figure')

    # Show the position of GC and anti-GC
    plt.plot(ra_gc_rad, dec_gc_rad, 'r+')
    plt.plot(ra_agc_rad, dec_agc_rad, 'r+')
    plt.text(ra_gc_rad, dec_gc_rad, "GC", color="r")
    # plt.text(ra_gc_rad, dec_gc_rad, "Anti-GC")

    if not fig_title is None:
        plt.title(fig_title, y=1.08)

    plt.grid(True)
    plt.subplots_adjust(top=0.95, bottom=0.05)

    if output is None:
        plt.show()
    else:
        plt.savefig(output)


# glide_plot([100, 0, 0])
# --------------------------------- END --------------------------------
