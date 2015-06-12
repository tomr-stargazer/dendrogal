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

from dendrogal.integrated_viewer import IntegratedViewer
from dendrogal.production.remove_degenerate_structures import reduce_catalog, selection_from_catalog
from dendrogal.hurt_image_overlay_demo import underlay_hurt_galaxy
from dendrogal.production.integrated_map_figure import IntegratedMapFigure

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'

# basic function: take in something analogous to an IntegratedViewer, and then let's figure out how to put like ellipses on it.

# this will require some WCS stuff  (WCSAxes almost certainly).

output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def make_quadrant_lbv_map(cloud_catalog, dendrogram, contour_select=True, 
                          alignment='lv', vscale=0.65, xscale=0.125, 
                          aspect=1/3, cmap='gray_r', **kwargs):

    # one: select things

    d = dendrogram
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

    lv_ells2 = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], 
                       width=2*cloud_catalog['major_sigma'][i]/l_scale_lbv, 
                       height=2*cloud_catalog['v_rms'][i]/v_scale_lbv) for i in range(len(cloud_catalog))]


    for e in lv_ells:
        iv.ax_lv.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor('black')
        e.set_linewidth(1)
        e.set_zorder(1)

    for e in lv_ells2:
        iv.ax_lv.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor('white')
        e.set_linewidth(2)
        e.set_zorder(0.95)

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

def make_quadrant_topdown_map(cloud_catalog, loc=None, figsize=(10,8)):

    fig = plt.figure(figsize=figsize)

    smalls = cloud_catalog[cloud_catalog['mass'] < 1e5]
    mediums = cloud_catalog[(cloud_catalog['mass'] < 1e6) & (cloud_catalog['mass'] > 1e5)]
    bigs = cloud_catalog[(cloud_catalog['mass'] < 1e7) & (cloud_catalog['mass'] > 1e6)]
    giants = cloud_catalog[cloud_catalog['mass'] > 1e7]

    plt.plot(smalls['y_sol'], smalls['x_sol'], 'co', ms=3, label=r"$10^4 M_\odot$")
    plt.plot(mediums['y_sol'], mediums['x_sol'], 'co', ms=5, label=r"$10^5 M_\odot$")
    plt.plot(bigs['y_sol'], bigs['x_sol'], 'co', ms=8, label=r"$10^6 M_\odot$")
    if len(giants) > 0:
        plt.plot(giants['y_sol'], giants['x_sol'], 'o', markerfacecolor='cyan', markeredgecolor='red', markeredgewidth=1.4, ms=11, label=r"$10^7 M_\odot$")

    leg = plt.legend(loc=loc, numpoints=1)
    leg.get_frame().set_alpha(0.5)

    underlay_hurt_galaxy(fig, u.kpc)

    return fig


def draw_solar_and_tangent_circles(fig):
    """ 
    Slightly hacky way to draw the two circles on the top-down map.

    """

    # assumes there's just one ax here

    ax = fig.get_axes()[0]

    R_0 = 8.34

    array_length = 100
    theta_array = np.linspace(0, 360, array_length)
    R_array_solar = np.repeat(R_0, array_length)
    R_array_tangent = np.repeat(R_0/2, array_length)

    x_solar_circle = R_array_solar * np.cos(np.radians(theta_array)) + R_0
    y_solar_circle = R_array_solar * np.sin(np.radians(theta_array))

    x_tangent_circle = R_array_tangent * np.cos(np.radians(theta_array)) + R_0/2
    y_tangent_circle = R_array_tangent * np.sin(np.radians(theta_array)) 

    lines1 = ax.plot(y_solar_circle, x_solar_circle, 'm', lw=2.5)

    lines2 = ax.plot(y_tangent_circle, x_tangent_circle, 'm--', lw=2)

    return lines1, lines2


def make_lv_map_new(cloud_catalog, dendrogram, ellipse_color=colorbrewer_red, 
                    ellipse_thickness_inner=0.9, ellipse_thickness_outer=1.9, 
                    **kwargs):

    d = dendrogram

    imf = IntegratedMapFigure(d.data, d.wcs, **kwargs)

    # A. DRAW CONTOURS - maybe not right now.
    # cloud_selection = selection_from_catalog(d, cloud_catalog)
    # mask = reduce(np.add, [structure.get_mask(subtree=True) for structure in cloud_selection])
    # imf.draw_contours(mask, colors=colorbrewer_red, linewidths=0.9, alpha=0.9)

    # B. DRAW ELLIPSES
    world_coordinates = np.vstack([cloud_catalog['x_cen'], cloud_catalog['y_cen'], cloud_catalog['v_cen']]).T
    lbv_pixels = d.wcs.wcs_world2pix(world_coordinates, 0)

    l_lbv_pixels = lbv_pixels[:,0]
    b_lbv_pixels = lbv_pixels[:,1]
    v_lbv_pixels = lbv_pixels[:,2]

    l_scale_lbv = imf.spatial_scale.value
    b_scale_lbv = imf.spatial_scale.value
    v_scale_lbv = imf.velocity_scale.value

    lv_ells = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], 
                       width=2*cloud_catalog['major_sigma'][i]/l_scale_lbv, 
                       height=2*cloud_catalog['v_rms'][i]/v_scale_lbv) for i in range(len(cloud_catalog))]

    lv_ells2 = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], 
                       width=2*cloud_catalog['major_sigma'][i]/l_scale_lbv, 
                       height=2*cloud_catalog['v_rms'][i]/v_scale_lbv) for i in range(len(cloud_catalog))]

    for e in lv_ells:
        imf.ax.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor(ellipse_color)
        e.set_linewidth(ellipse_thickness_inner)
        e.set_zorder(0.95)

    for e in lv_ells2:
        imf.ax.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor('white')
        e.set_linewidth(ellipse_thickness_outer)
        # e.set_alpha(0.8)
        e.set_zorder(0.9)

    lon = imf.ax.coords['glon']
    lon.set_ticks(spacing=5*u.deg, color='white', exclude_overlapping=True)
    lon.display_minor_ticks(True)
    lon.set_axislabel(r"$l$ (deg)", minpad=1.5)

    vlsr = imf.ax.coords['vopt']
    vlsr.set_ticks(spacing=25*u.m/u.s, color='white', exclude_overlapping=True) # erroneous units - why!?
    vlsr.display_minor_ticks(True)
    vlsr.set_axislabel(r"$v_{LSR}$ (km s$^{-1}$)")
    vlsr.set_ticklabel_position('lr')


    imf.fig.canvas.draw()

    return imf