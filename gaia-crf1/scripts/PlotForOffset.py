# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 09:20:26 2016

@author: Neo
"""

import matplotlib.pyplot as plt

## plot 
### d_RA - RA
#plt.errorbar(RAg/3600, d_RA, yerr=e_dRA, fmt='.')
#plt.xlim([0, 360])
#plt.xticks(np.arange(0, 400, 90))
#plt.xlabel('$R.A.(^o)$')
#plt.ylabel('$\Delta R.A.^*(mas)$')
#plt.axhline(xmax=360, linewidth=1, color='r')
#plt.savefig('dRA_RA.eps', dpi=100)
#plt.close()
### d_RA - DE
#plt.errorbar(DEg/3600, d_RA, yerr=e_dRA, fmt='.')
#plt.xlim([-90, 90])
#plt.xticks(np.arange(-90, 100, 30))
#plt.xlabel('$Dec.(^o)$')
#plt.ylabel('$\Delta R.A.^*(mas)$')
#plt.axhline(xmin=-90, xmax=90, linewidth=1, color='r')
#plt.savefig('dRA_DE.eps', dpi=100)
#plt.close()
### d_DE - RA
#plt.errorbar(RAg/3600, d_DE, yerr=e_dDE, fmt='.')
#plt.xlim([0, 360])
#plt.xticks(np.arange(0, 400, 90))
#plt.xlabel('$R.A.(^o)$')
#plt.ylabel('$\Delta Dec.(mas)$')
#plt.axhline(xmax=360, linewidth=1, color='r')
#plt.savefig('dDE_RA.eps', dpi=100)
#plt.close()
### d_DE - DE
#plt.errorbar(DEg/3600, d_DE, yerr=e_dDE, fmt='.')
#plt.xlim([-90, 90])
#plt.xticks(np.arange(-90, 100, 30))
#plt.xlabel('$Dec.(^o)$')
#plt.ylabel('$\Delta Dec.(mas)$')
#plt.axhline(xmin=-90, xmax=90, linewidth=1, color='r')
#plt.savefig('dDE_DE.eps', dpi=100)
#plt.close()
### e_RAi - e_RAg
#plt.plot(e_RAi, e_RAg, '.')
#plt.xlabel('eRA_icrf2(mas)')
#plt.ylabel('eRA_gaia(mas)')
#plt.xlim([0, 20])
#plt.ylim([0, 20])
##plt.show()
#plt.savefig('eRAi_eRAg.eps', dpi=100)
#plt.close()
### e_DEi - e_DEg
#plt.plot(e_DEi, e_DEg, '.')
#plt.xlabel('eDE_icrf2(mas)')
#plt.ylabel('eDE_gaia(mas)')
#plt.xlim([0, 20])
#plt.ylim([0, 20])
##plt.show()
#plt.savefig('eDEi_eDEg.eps', dpi=100)
#plt.close()
#
## A complex plot 
## scatter of dRA, dDE
#x = d_RA
#y = d_DE
#from matplotlib.ticker import NullFormatter
#
#nullfmt = NullFormatter()         # no labels
#
## definitions for the axes
#left, width = 0.1, 0.65
#bottom, height = 0.1, 0.65
#bottom_h = left_h = left + width + 0.02
#
#rect_scatter = [left, bottom, width, height]
#rect_histx = [left, bottom_h, width, 0.2]
#rect_histy = [left_h, bottom, 0.2, height]
#
## start with a rectangular Figure
#plt.figure(1, figsize=(8, 8))
#
#axScatter = plt.axes(rect_scatter)
#axHistx = plt.axes(rect_histx)
#axHisty = plt.axes(rect_histy)
#
## no labels
#axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)
#
## the scatter plot:
#axScatter.scatter(x, y)
#axScatter.set_xlabel('$\Delta R.A.^*(mas)$')
#axScatter.set_ylabel('$\Delta Dec.(mas)$')
## now determine nice limits by hand:
#binwidth = 1.0
#xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
##lim = (int(xymax/binwidth) + 1) * binwidth
#lim = 10.0
#
#axScatter.set_xlim((-lim, lim))
#axScatter.set_ylim((-lim, lim))
#
#bins = np.arange(-lim, lim + binwidth, binwidth)
#axHistx.hist(x, bins=bins)
#axHisty.hist(y, bins=bins, orientation='horizontal')
#axHistx.set_xlim(axScatter.get_xlim())
#axHisty.set_ylim(axScatter.get_ylim())
#axHistx.set_ylim([0,800])
#axHistx.set_yticks(range(0, 800, 200))
#axHisty.set_xlim([0,800])
#axHisty.set_xticks(range(200, 800, 200))
#plt.savefig('dRA_dDE.eps', dpi=100)
#plt.close()

## plot the positional difference.
#scale = np.sqrt(20)
#for i in range(len(ICRFg)):
#    d = np.sqrt(d_RA[i]**2 + d_DE[i]**2)
#    if d > 10:
#        plt.arrow(RAg[i]/3600, DEg[i]/3600, \
#        d_RA[i]/d*10*scale, d_DE[i]/d*10*scale, color='r')
#    else:
#        plt.arrow(RAg[i]/3600, DEg[i]/3600, \
#        d_RA[i]*scale, d_DE[i]*scale, color='k')
##    if ICRFg[i] in icrf2d:
##        plt.plot(RAg[i]/3600, DEg[i]/3600, 'b.')
##    else:
##        plt.plot(RAg[i]/3600, DEg[i]/3600, 'r.')
#    plt.plot(RAg[i]/3600, DEg[i]/3600, 'b.')
#        
#plt.axhline(y=-80, xmin=10, xmax=30)
#plt.text(20,-60, '1mas')
#plt.xlim([0, 360])
#plt.ylim([-90, 90])
#plt.xticks(np.arange(0, 400, 30))
#plt.yticks(np.arange(-90, 100, 30))
#plt.show()