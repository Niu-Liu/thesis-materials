#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: glide_calc.py
"""
Created on Mon Mar 26 16:21:26 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Some functions related to glide component calculation.

History
N.Liu, 20 Mar 2018:

"""

import numpy as np
from numpy import cos, sin


__all__ = ["glide_gen", "glide_calc", "GA_glide_decomposed",
           "glide_field_gen"]


# -----------------------------  FUNCTIONS -----------------------------
def glide_gen(g, RAdeg, DCdeg, err=None):
    '''Given apex and amplitude, calculate the glide vector.

    NOTE: inputs RAdeg/DCdeg should be given in degree.

    Parameters
    ----------
    g : float
        amplitude of glide vector, default value is 1.0.
    RAdeg/DCdeg : float
        direction of apex expressed in right ascension and declination
        given in degree. For Galactic center, they could be set as
        267.0/-29.0.
    err : array of float
        uncertainties of g, RAdeg and DCdeg given in rad/deg/deg,
        default value is None.

    Returns
    ----------
    [g1, g2, g3] : array of float
        glide vector
    errn : array of float
        vector of glide formal uncertainties
    '''

    # Degree -> radian
    RArad = np.deg2rad(RAdeg)
    DCrad = np.deg2rad(DCdeg)

    g1 = g * cos(RArad) * cos(DCrad)
    g2 = g * sin(RArad) * cos(DCrad)
    g3 = g * sin(DCrad)

    if err is None:
        return np.array([g1, g2, g3])
    else:
        # Calculate the uncertainties.
        M = np.array([
            [cos(RArad) * cos(DCrad),
             -g * sin(RArad) * cos(DCrad),
             -g * cos(RArad) * sin(DCrad)],
            [sin(RArad) * cos(DCrad),
             g * cos(RArad) * cos(DCrad),
             -g * sin(RArad) * sin(DCrad)],
            [sin(DCrad),
             0,
             g * cos(DCrad)]])

        # degree -> radian
        err[1] = np.deg2rad(err[1])
        err[2] = np.deg2rad(err[2])

        errn = np.sqrt(np.dot(M**2, err**2))

        return np.array([g1, g2, g3]), errn


def glide_calc(gv, err_gv=None):
    '''Calculate the apex and amplitude for a given glide vector.

    Parameters
    ----------
    gv : array of float
        glide vector
    err_gv : array of float
        formal error of glide vector, default value is None

    Returns
    ----------
    g/err_g : amplitude and its uncertainty
    RAdeg/err_RA : right ascension and its uncertainty in degree
    DCdeg/err_DC : declination and its uncertainty in degree
    '''

    # Calculate the parameters.
    # Amplitude
    g = np.sqrt(np.sum(gv ** 2))

    # Right ascension
    RArad = np.arctan2(gv[1], gv[0])
    # arctan2 function result is always between -pi and pi.
    # However, RA should be between 0 and 2*pi.
    if RArad < 0:
        RArad += 2 * np.pi

    # Declination
    DCrad = np.arctan2(gv[2], np.sqrt(gv[0]**2 + gv[1]**2))

    # radian -> degree
    RAdeg = np.rad2deg(RArad)
    DCdeg = np.rad2deg(DCrad)

    if err_gv is None:
        return g, RAdeg, DCdeg
    else:
        # # Method 1)
        # M = np.array([
        #     [cos(RArad) * cos(DCrad),
        #      g * sin(RArad) * cos(DCrad),
        #      g * cos(RArad) * sin(DCrad)],
        #     [sin(RArad) * cos(DCrad),
        #      g * cos(RArad) * cos(DCrad),
        #      g * sin(RArad) * sin(DCrad)],
        #     [sin(DCrad),
        #      0,
        #      g * cos(DCrad)]])
        # M_1 = np.linalg.inv(M**2)
        # errn = np.sqrt(np.dot(M_1, err_gv**2))
        # errA, errRAr, errDCr = errn
        # However, the obtained error are always meaningless (minus
        # value inner the square operation).

        # #) Method 2)
        par_g = gv / g
        gxy2 = gv[0]**2 + gv[1]**2
        gxy = np.sqrt(gxy2)
        par_RA = np.array([-gv[1], gv[0], 0]) / gxy2
        par_DC = np.array([-gv[0] * gv[2],
                           -gv[1] * gv[2],
                           gxy2]) / g**2 / gxy
        errA = np.sqrt(np.dot(par_g**2, err_gv**2))
        errRAr = np.sqrt(np.dot(par_RA**2, err_gv**2))
        errDCr = np.sqrt(np.dot(par_DC**2, err_gv**2))

        # --------------- 20 Mar 2018 ----------------------------------
        # errA = np.sqrt(np.dot(err_gv, err_gv))
        # --------------- 20 Mar 2018 ----------------------------------

        # radian -> degree
        errRA = np.rad2deg(errRAr)
        errDC = np.rad2deg(errDCr)

        return g, RAdeg, DCdeg, errA, errRA, errDC


def GA_glide_decomposed(gv, err_gv):
    """Given a glide vector, get the GA component and non-GA component.


    'G' stands for amplitude while 'g' for vector

    Parameters
    ----------
    gv : array of float
        glide vector
    err_gv : array of float
        formal error of glide vector

    Returns
    ----------
    G_GA/err_GA : amplitude and its uncertainty of GA component
    g_NonGA : non-GA glide component
    RA_NonGA/DC_NonGA : apex of non-GA component
    """

    GA_hat = glide_gen(1.0, 266.4, -28.9)

    # GA component
    G_GA = np.dot(gv, GA_hat)
    errG_GA = np.dot(err_gv, GA_hat)

    # non-GA component
    g_NonGA = gv - G_GA * GA_hat
    err_NonGA = err_gv - errG_GA * GA_hat
    [G_NonGA, RA_NonGA, DC_NonGA,
     errG_NonGA, errRA_NonGA, errDC_NonGA] = glide_calc(g_NonGA, err_NonGA)

    return [G_GA, errG_GA, G_NonGA, RA_NonGA, DC_NonGA,
            errG_NonGA, errRA_NonGA, errDC_NonGA]


def test_code(flag):
    '''Used for verified the code.
    '''

    if flag is 1:
        # print(glide_calc(np.array([1.2, 2.4, 7.5]),
                         # np.array([0.6, 0.8, 1.0])))

        # Verify the result of Table 1 in Titov et al.(2011) (A&A 529, A91)
        print('\nTo verify the result of Table 1 in Titov et al.(2011):')
        # Result in the paper
        dip = np.array([[-0.7, -5.9, -2.2],
                        [-1.9, -7.1, -3.5],
                        [-1.6, -6.0, -2.9],
                        [-3.4, -6.2, -3.8],
                        [+4.2, -4.6, +1.4]])
        derr = np.array([[0.8, 0.9, 1.0],
                         [1.0, 1.1, 1.2],
                         [0.9, 1.0, 1.1],
                         [1.0, 1.2, 1.2],
                         [1.3, 1.3, 1.7]])
        par = np.array([[6.4, 263, -20],
                        [8.1, 255, -25],
                        [6.9, 255, -25],
                        [8.0, 241, -28],
                        [6.4, 312, 13]])
        perr = np.array([[1.5, 11, 12],
                         [1.9, 11, 11],
                         [1.7, 11, 11],
                         [2.0, 12, 12],
                         [2.5, 16, 21]])
        for (dipi, derri, pari, perri) in zip(dip, derr, par, perr):
            print('Result in paper: ', pari, perri)
            print('Result of my calculation:', glide_calc(dipi, derri))
            print('For g:', np.sqrt(np.dot(derri, derri)))

        # Verify the result of Table 1 in Liu et al.(2012) (A&A 548, A50)
        print('\nTo verify the result of Table 1 in Liu et al.(2012):')
        # Result in the paper
        dip = np.array([[+1.07, -0.20, +0.18],
                        [+0.22, -0.03, +0.04],
                        [+0.01, -0.14, +0.28],
                        [+0.89, -0.08, +0.04]])
        derr = np.array([[0.08, 0.08, 0.09],
                         [0.05, 0.05, 0.05],
                         [0.07, 0.07, 0.07],
                         [0.02, 0.02, 0.02]])
        par = np.array([[1.10],
                        [0.22],
                        [0.31],
                        [0.89]])
        perr = np.array([[0.14],
                         [0.09],
                         [0.12],
                         [0.03]])
        for (dipi, derri, pari, perri) in zip(dip, derr, par, perr):
            print('Result in paper: ', pari, perri)
            print('Result of my calculation:', glide_calc(dipi, derri))
            print('For g:', np.sqrt(np.dot(derri, derri)))

        # Verify the result of Table 1 in Titov et al.(2013) (A&A 559, A95)
        print('\nTo verify the result of Table 1 in Titov et al.(2013):')
        # Result in the paper
        dip = np.array([[-0.4, -5.7, -2.8],
                        [+0.7, -6.2, -3.3],
                        [+0.7, -6.2, -3.3]])
        derr = np.array([[0.7, 0.8, 0.9],
                         [0.8, 0.9, 1.0],
                         [0.9, 1.0, 1.0]])
        par = np.array([[6.4, 266, -26],
                        [7.1, 277, -28],
                        [7.1, 277, -28]])
        perr = np.array([[1.1, 7, 7],
                         [1.3, 7, 7],
                         [1.4, 9, 8]])
        for (dipi, derri, pari, perri) in zip(dip, derr, par, perr):
            print('Result in paper: ', pari, perri)
            print('Result of my calculation:', glide_calc(dipi, derri))
            print('For g:', np.sqrt(np.dot(derri, derri)))

        # Verify the result of Table 1 in Titov et al.(2018) (A&A 610, A36)
        print('\nTo verify the result of Table 1 in Titov et al.(2018):')
        # Result in the paper
        dip = np.array([[-0.7, -5.9, -2.2],
                        [-2.6, -5.1, -1.1],
                        [-0.4, -5.7, -2.8],
                        [-0.3, -5.4, -1.1],
                        [+0.2, -3.3, -4.9]])
        derr = np.array([[0.8, 0.9, 1.0],
                         [0.3, 0.3, 0.4],
                         [0.7, 0.8, 0.9],
                         [0.3, 0.3, 0.5],
                         [0.8, 0.8, 1.0]])
        par = np.array([[6.4, 263, -20],
                        [5.8, 243, -11],
                        [6.4, 266, -26],
                        [5.6, 267, -11],
                        [5.9, 273, -56]])
        perr = np.array([[1.5, 11, 12],
                         [0.4, 4, 4],
                         [1.1, 7, 7],
                         [0.4, 3, 3],
                         [1.0, 13, 9]])
        for (dipi, derri, pari, perri) in zip(dip, derr, par, perr):
            print('Result in paper: ', pari, perri)
            print('Result of my calculation:', glide_calc(dipi, derri))
            print('For g:', np.sqrt(np.dot(derri, derri)))

        # Verify the result in David Mayer(draft)
        print('\nTo verify the result of Table 1 in draft of David Mayer:')
        # Result in the paper
        dip = np.array([[+2., +26., +32.],
                        [+5., +25., +44.],
                        [+1., +2., +21.],
                        [+3., -8., +27.],
                        [-2., +5., +6.],
                        [-3., -8., +19.]])
        derr = np.array([[15., 16., 14.],
                         [15., 16., 14.],
                         [15., 16., 14.],
                         [15., 16., 14.],
                         [15., 16., 13.],
                         [15., 16., 14.]])
        par = np.array([[42, 85, +51],
                        [51, 78, +60],
                        [21, 77, +83],
                        [29, 291, +73],
                        [8, 108, +53],
                        [21, 250, +65]])
        perr = np.array([[15, 33, 19],
                         [15, 35, 16],
                         [14, 350, 41],
                         [14, 102, 32],
                         [14, 175, 107],
                         [13, 104, 43]])
        for (dipi, derri, pari, perri) in zip(dip, derr, par, perr):
            print('Result in paper: ', pari, perri)
            print('Result of my calculation:', glide_calc(dipi, derri))
            print('For g:', np.sqrt(np.dot(derri, derri)))

    elif flag is 2:
        g, RA, DC = 5.4, 256.0, -30
        err_g, err_RA, err_DC = 0.6, 4.0, 5.0
        print('Input: ', g, RA, DC, err_g, err_RA, err_DC)
        print('Radian for angles', np.deg2rad(err_RA), np.deg2rad(err_DC))
        gv, err_gv = glide_gen(g, RA, DC,
                               np.array([err_g, err_RA, err_DC]))
        g1, RA1, DC1, err_g1, err_RA1, err_DC1 = glide_calc(gv, err_gv)
        print('Output: ', g1, RA1, DC1, err_g1, err_RA1, err_DC1)


# --------------------------------- MAIN --------------------------------
if __name__ == '__main__':
    test_code(1)
# --------------------------------- END --------------------------------
