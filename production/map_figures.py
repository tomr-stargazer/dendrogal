""" This is a script that will help us create figure maps. """

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from ..integrated_viewer import IntegratedViewer
from .remove_degenerate_structures import reduce_catalog, selection_from_catalog
from ..hurt_image_overlay_demo import underlay_hurt_galaxy

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'

# basic function: take in something analogous to an IntegratedViewer, and then let's figure out how to put like ellipses on it.

# this will require some WCS stuff  (WCSAxes almost certainly).

output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def make_quadrant_lbv_map(cloud_catalog, dendgrogram, contour_select=True, 
                          alignment='lv', vscale=0.65, xscale=0.125, 
                          aspect=1/3, cmap='gray_r', **kwargs):

    # one: select things

    d = dendgrogram
    dv = d.viewer()

    cloud_selection = selection_from_catalog(d, cloud_catalog)

    if contour_select:
        dv.hub.select(2, cloud_selection, subtree=False)
        dv.hub.select(1, None)
        dv.hub.select(3, None)

    iv = IntegratedViewer(d, dv.hub, cmap=cmap, clip_velocity=False, 
                          alignment=alignment, aspect=aspect, **kwargs)

    # let's make ellipses
    # copied largely from comparison_to_other_catalogs.py

    world_coordinates = np.vstack([cloud_catalog['x_cen'], cloud_catalog['y_cen'], cloud_catalog['v_cen']]).T

    lbv_pixels = iv.dendrogram.wcs.wcs_world2pix(world_coordinates, 0)

    l_lbv_pixels = lbv_pixels[:,0]
    b_lbv_pixels = lbv_pixels[:,1]
    v_lbv_pixels = lbv_pixels[:,2]

    # l_scale_lv = (l_lv_pixels[1] - l_lv_pixels[0])/(lcen_column[1]-lcen_column[0])
    # v_scale_lv = (v_lv_pixels[1] - v_lv_pixels[0])/(vcen_column[1]-vcen_column[0])
    l_scale_lbv = xscale # 0.125 in quadrant 1
    b_scale_lbv = xscale # 0.125 in quadrant 1
    v_scale_lbv = vscale # 0.65 in quadrant 1

    lv_ells = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], 
                       width=2*cloud_catalog['major_sigma'][i]/l_scale_lbv, 
                       height=2*cloud_catalog['v_rms'][i]/v_scale_lbv) for i in range(len(cloud_catalog))]
    for e in lv_ells:
        iv.ax_lv.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor(colorbrewer_red)
        e.set_linewidth(1.5)
        e.set_zorder(1)

    lon = iv.ax_lv.coords['glon']
    lon.set_ticks(spacing=5*u.deg)
    lon.display_minor_ticks(True)
    lon.set_axislabel(r"$l$ (deg)", minpad=1.5)

    vlsr = iv.ax_lv.coords['vopt']
    vlsr.set_ticks(spacing=25*u.m/u.s) # erroneous units - why!?
    vlsr.display_minor_ticks(True)

    vlsr.set_axislabel(r"$v_{LSR}$ (km s$^{-1}$)")            

    if alignment != 'lv':
        lb_ells = [(Ellipse(xy=zip(l_lbv_pixels, b_lbv_pixels)[i], 
                           angle=cloud_catalog['position_angle'][i],
                           width=2*cloud_catalog['major_sigma'][i]/l_scale_lbv, 
                           height=2*cloud_catalog['minor_sigma'][i]/b_scale_lbv),
                    cloud_catalog['distance'][i]) for i in range(len(cloud_catalog))]

        for e, distance in lb_ells:
            iv.ax_lb.add_artist(e)
            e.set_facecolor(colorbrewer_red)
            e.set_alpha(0.3)
            e.set_edgecolor(colorbrewer_red)
            e.set_linewidth(1.5)
            e.set_zorder(2 - 0.001 * distance)

        iv.ax_lb.set_xlabel(r"Galactic longitude $l$ (deg)")
        iv.ax_lb.set_ylabel(r"Galactic latitude $b$ (deg)")

    iv.fig.canvas.draw()

    return iv

def make_quadrant_topdown_map(cloud_catalog, loc=None):

    fig = plt.figure()

    smalls = cloud_catalog[cloud_catalog['mass'] < 1e5]
    mediums = cloud_catalog[(cloud_catalog['mass'] < 1e6) & (cloud_catalog['mass'] > 1e5)]
    bigs = cloud_catalog[(cloud_catalog['mass'] < 1e7) & (cloud_catalog['mass'] > 1e6)]
    giants = cloud_catalog[cloud_catalog['mass'] > 1e7]

    plt.plot(smalls['y_sol'], smalls['x_sol'], 'co', ms=4, label=r"$10^4 M_\odot$")
    plt.plot(mediums['y_sol'], mediums['x_sol'], 'co', ms=7, label=r"$10^5 M_\odot$")
    plt.plot(bigs['y_sol'], bigs['x_sol'], 'co', ms=9, label=r"$10^6 M_\odot$")
    if len(giants) > 0:
        plt.plot(giants['y_sol'], giants['x_sol'], 'co', ms=10, label=r"$10^7 M_\odot$")

    leg = plt.legend(loc=loc, numpoints=1)
    leg.get_frame().set_alpha(0.5)

    underlay_hurt_galaxy(fig, u.kpc)

    return fig

