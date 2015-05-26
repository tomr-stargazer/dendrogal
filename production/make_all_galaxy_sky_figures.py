"""
Makes the all-Galaxy sky figure.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt

path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def make_that_figure(catalog):

    fig = plt.figure(figsize=(8,5))

    ax_lb = fig.add_subplot(611)
    ax_lv = fig.add_subplot(312, sharex=ax_lb)
    ax_ld = fig.add_subplot(313, sharex=ax_lb)

    ax_lb.plot(catalog['x_cen'], catalog['y_cen'], 'k.', ms=2)
    ax_lb.plot(catalog['x_cen']-360, catalog['y_cen'], 'k.', ms=2)

    ax_lv.plot(catalog['x_cen'], catalog['v_cen'], 'k.', ms=2)    
    ax_lv.plot(catalog['x_cen']-360, catalog['v_cen'], 'k.', ms=2)    

    ax_ld.errorbar(
        catalog['x_cen'], catalog['distance'], 
        yerr=[catalog['error_distance_minus'], catalog['error_distance_plus']], 
        fmt='k,', capsize=0)

    ax_ld.errorbar(
        catalog['x_cen']-360, catalog['distance'], 
        yerr=[catalog['error_distance_minus'], catalog['error_distance_plus']], 
        fmt='k,', capsize=0)

    ax_lb.set_xlim(180, -180)
    ax_lb.set_aspect('equal')
    ax_lb.set_ylabel("$b$", rotation='horizontal')
    ax_lb.set_xlabel("$l$")

    ax_lb.set_ylim(-10,10)

    ax_lv.set_ylim(-200, 200)
    ax_lv.set_ylabel("$v_{LSR}$ (km/s)")

    x_ticks = np.arange(-180, 180+30, 30)

    x_tick_labels = [x if x>=0 else x+360 for x in x_ticks]

    ax_lb.set_xticks(x_ticks)
    ax_lb.set_xticklabels(x_tick_labels)

    ax_ld.set_ylabel("Distance from Sun")
    ax_ld.set_xlabel("$l$")

    return fig


def save_that_figure(catalog):

    make_that_figure(catalog).savefig(path+'all_galaxy_longitude_trio.pdf')
