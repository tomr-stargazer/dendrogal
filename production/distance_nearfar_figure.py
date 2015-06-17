"""
Makes the figure with the distance near/far stuff.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt

output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def make_nearfar_figure(catalog):

    fig = plt.figure(figsize=(13,4))

    ax1 = fig.add_subplot(131)

    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    far_catalog = catalog[ catalog['KDA_resolution'] =='F']

    near_catalog = catalog[ catalog['KDA_resolution'] =='N']

    ambig_catalog = catalog[ catalog['KDA_resolution'] == 'A']

    ax1.plot(catalog['p_far'], catalog['p_near'], 'k.')    
    ax1.plot(near_catalog['p_far'], near_catalog['p_near'], 'b.')
    ax1.plot(far_catalog['p_far'], far_catalog['p_near'], 'r.')

    ax1.set_xlabel("p_far", fontsize=16)
    ax1.set_ylabel("p_near", fontsize=16)

    ax2.plot(np.log10(far_catalog['far_size']), np.log10(far_catalog['v_rms']), 'r.')
    ax2.plot(np.log10(far_catalog['near_size']), np.log10(far_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='r', mew=0.9, ms=6)

    ax2.plot(np.log10(near_catalog['far_size']), np.log10(near_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='b', mew=0.9, ms=6)
    ax2.plot(np.log10(near_catalog['near_size']), np.log10(near_catalog['v_rms']), 'b.')

    ax2.set_xlabel("log(R/pc)", fontsize=16)
    ax2.set_ylabel("log($\sigma_v$)", fontsize=16)

    ax3.plot(near_catalog['near_distance'], near_catalog['near_z_gal']*1000, 'b.')
    ax3.plot(near_catalog['far_distance'], near_catalog['far_z_gal']*1000, 'x', markerfacecolor='none', markeredgecolor='b', mew=0.9, ms=6)

    ax3.plot(far_catalog['near_distance'], far_catalog['near_z_gal']*1000, 'x', markerfacecolor='none', markeredgecolor='r', mew=0.9, ms=6)
    ax3.plot(far_catalog['far_distance'], far_catalog['far_z_gal']*1000, 'r.')

    ax3.set_xlabel("$D_\odot$ (kpc)", fontsize=16)
    ax3.set_ylabel("$z_{gal}$ (pc)", fontsize=16)

    # ax2.plot(np.log10(ambig_catalog['far_size']), np.log10(ambig_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='k', mew=0.5, ms=3.5)
    # ax2.plot(np.log10(ambig_catalog['near_size']), np.log10(ambig_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='k', mew=0.5, ms=3.5)

    fig.canvas.draw()
    fig.tight_layout(pad=0.5)
    fig.canvas.draw() 

    return fig

def save_nearfar_figure(catalog):
    fig = make_nearfar_figure(catalog)
    fig.savefig(output_path+"nearfar_figure.pdf")