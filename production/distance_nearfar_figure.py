"""
Makes the figure with the distance near/far stuff.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt

output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")


def make_nearfar_figure(catalog):

    fig = plt.figure(figsize=(7,9))

    ax1 = fig.add_subplot(321)

    ax2 = fig.add_subplot(323)
    ax2b = fig.add_subplot(324)

    ax3 = fig.add_subplot(325)
    ax3b = fig.add_subplot(326)

    far_catalog = catalog[ catalog['KDA_resolution'] =='F']

    near_catalog = catalog[ catalog['KDA_resolution'] =='N']

    ambig_catalog = catalog[ catalog['KDA_resolution'] == 'A']

    ax1.plot(catalog['p_far'], catalog['p_near'], 'k.')    
    ax1.plot(near_catalog['p_far'], near_catalog['p_near'], 'b.')
    ax1.plot(far_catalog['p_far'], far_catalog['p_near'], 'r.')

    ax1.set_xlabel("p_far", fontsize=14)
    ax1.set_ylabel("p_near", fontsize=14)

    ###############

    ax2.plot(np.log10(far_catalog['far_size']), np.log10(far_catalog['v_rms']), 'r.', zorder=1, ms=6)
    ax2.plot(np.log10(far_catalog['near_size']), np.log10(far_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='0.75', mew=0.9, ms=4, zorder=0.5)

    ax2.plot(np.log10(near_catalog['far_size']), np.log10(near_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='0.75', mew=0.9, ms=4, zorder=0.5)
    ax2.plot(np.log10(near_catalog['near_size']), np.log10(near_catalog['v_rms']), 'b.', zorder=1, ms=6)

    ax2.set_xlabel("log(R/pc)", fontsize=14)
    ax2.set_ylabel("log($\sigma_v$)", fontsize=14)

    ax2.set_ylim(0, 1.4)

    ###############

    ax2b.plot(np.log10(far_catalog['far_size']), np.log10(far_catalog['v_rms']), '.', color='0.75', zorder=0.5)
    ax2b.plot(np.log10(far_catalog['near_size']), np.log10(far_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='r', mew=0.9, ms=4, zorder=1)

    ax2b.plot(np.log10(near_catalog['far_size']), np.log10(near_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='b', mew=0.9, ms=4, zorder=1)
    ax2b.plot(np.log10(near_catalog['near_size']), np.log10(near_catalog['v_rms']), '.', color='0.75', zorder=0.5)

    ax2b.set_xlabel("log(R/pc)", fontsize=14)
    # ax2b.set_ylabel("log($\sigma_v$)", fontsize=14)

    ax2b.set_ylim(0, 1.4)

    ###############

    ax3.plot(near_catalog['near_distance'], near_catalog['near_z_gal']*1000, 'b.', zorder=1)
    ax3.plot(near_catalog['far_distance'], near_catalog['far_z_gal']*1000, 'x', markerfacecolor='none', markeredgecolor='0.75', mew=0.9, ms=4, zorder=0.5)

    ax3.plot(far_catalog['near_distance'], far_catalog['near_z_gal']*1000, 'x', markerfacecolor='none', markeredgecolor='0.75', mew=0.9, ms=4, zorder=0.5)
    ax3.plot(far_catalog['far_distance'], far_catalog['far_z_gal']*1000, 'r.', zorder=1)

    ax3.set_xlabel("$d_\odot$ (kpc)", fontsize=14)
    ax3.set_ylabel("$z_{gal}$ (pc)", fontsize=14)

    ###############

    ax3b.plot(near_catalog['near_distance'], near_catalog['near_z_gal']*1000, '.', color='0.75', zorder=0.5)
    ax3b.plot(near_catalog['far_distance'], near_catalog['far_z_gal']*1000, 'x', markerfacecolor='none', markeredgecolor='b', mew=0.9, ms=4, zorder=1)

    ax3b.plot(far_catalog['near_distance'], far_catalog['near_z_gal']*1000, 'x', markerfacecolor='none', markeredgecolor='r', mew=0.9, ms=4, zorder=1)
    ax3b.plot(far_catalog['far_distance'], far_catalog['far_z_gal']*1000, '.', color='0.75', zorder=0.5)

    ax3b.set_xlabel("$d_\odot$ (kpc)", fontsize=14)
    # ax3b.set_ylabel("$z_{gal}$ (pc)", fontsize=14)

    # ax2.plot(np.log10(ambig_catalog['far_size']), np.log10(ambig_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='k', mew=0.5, ms=3.5)
    # ax2.plot(np.log10(ambig_catalog['near_size']), np.log10(ambig_catalog['v_rms']), 'x', markerfacecolor='none', markeredgecolor='k', mew=0.5, ms=3.5)

    ax1.text(0.05, 0.91, '(a)', transform=ax1.transAxes, fontsize=14, family='serif')
    ax2.text(0.05, 0.91, '(b)', transform=ax2.transAxes, fontsize=14, family='serif')
    ax2b.text(0.05, 0.91, '(c)', transform=ax2b.transAxes, fontsize=14, family='serif')
    ax3.text(0.05, 0.91, '(d)', transform=ax3.transAxes, fontsize=14, family='serif')
    ax3b.text(0.05, 0.91, '(e)', transform=ax3b.transAxes, fontsize=14, family='serif')

    ax2.text(0.05, 0.825, 'assigned distance', transform=ax2.transAxes, fontsize=11, family='serif')
    ax2b.text(0.05, 0.825, 'rejected distance', transform=ax2b.transAxes, fontsize=11, family='serif')

    ax3.text(0.05, 0.825, 'assigned distance', transform=ax3.transAxes, fontsize=11, family='serif')
    ax3b.text(0.05, 0.825, 'rejected distance', transform=ax3b.transAxes, fontsize=11, family='serif')

    fig.canvas.draw()
    fig.tight_layout(pad=0.65)
    fig.canvas.draw() 

    return fig

def save_nearfar_figure(catalog):
    fig = make_nearfar_figure(catalog)
    fig.savefig(output_path+"nearfar_figure.pdf", bbox_inches='tight')