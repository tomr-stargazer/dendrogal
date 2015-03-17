"""
A script allowing (visual, statistical) comparisons between our catalog and others.

Some included here:
Dame+1986 (First Quadrant, 33 objects)
Solomon+1987 (First Quadrant, 464 objects)
Scoville+1987 (First Quadrant, 1427 objects)


"""

from __future__ import division

import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from astropy import table

data_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/other_catalogs/")

dame86_catalog = table.Table.read(data_path+'Dame+86_catalog_sanitized.txt', format='ascii', header_start=2)
garcia_catalog = table.Table.read(data_path+'Garcia_catalog.fit')

id_column = dame86_catalog['ID']
vcen_column = np.array([x.split(',')[1] for x in id_column]).astype('float')
lcen_column = dame86_catalog['LII']
sky_radius_column = np.degrees(dame86_catalog['R']/(dame86_catalog['D']*1000))

garcia_sky_radius = np.degrees(garcia_catalog['R']/(garcia_catalog['D']*1000))

def make_figures_to_show_lv_stuff():
    fig1 = plt.figure()

    plt.plot(lcen_column, vcen_column, 'bo')
    plt.gca().invert_xaxis()

    plt.show()

    fig2 = plt.figure()

    ax = fig2.add_subplot(111)

    ells = [Ellipse(xy=zip(lcen_column, vcen_column)[i], 
                    width=2*sky_radius_column[i], height=dame86_catalog['DV'][i]) for i in range(len(dame86_catalog))]
    for e in ells:
        ax.add_artist(e)
        e.set_facecolor('none')
    plt.xlim(65,0)
    plt.ylim(0,150)
    plt.show()

    return fig1, fig2

def plot_dame_ellipses_on_integrated_viewer(integrated_viewer):

    iv = integrated_viewer

    world_coordinates = np.vstack([lcen_column, np.repeat(0, len(lcen_column)), vcen_column]).T

    lv_pixels = iv.dendrogram.wcs.wcs_world2pix(world_coordinates, 0)

    l_lv_pixels = lv_pixels[:,0]
    v_lv_pixels = lv_pixels[:,2]

    # l_scale_lv = (l_lv_pixels[1] - l_lv_pixels[0])/(lcen_column[1]-lcen_column[0])
    # v_scale_lv = (v_lv_pixels[1] - v_lv_pixels[0])/(vcen_column[1]-vcen_column[0])
    l_scale_lv = 0.125
    v_scale_lv = 0.65

    lv_ells = [Ellipse(xy=zip(l_lv_pixels, v_lv_pixels)[i], 
                       width=2*sky_radius_column[i]/l_scale_lv, 
                       height=dame86_catalog['DV'][i]/v_scale_lv) for i in range(len(dame86_catalog))]
    for e in lv_ells:
        iv.ax_lv.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor('blue')


def plot_garcia_ellipses_on_iv(integrated_viewer):

    iv = integrated_viewer

    world_coordinates = np.vstack([garcia_catalog['GLON'].data, garcia_catalog['GLAT'].data, garcia_catalog['Vlsr'].data]).T

    lv_pixels = iv.dendrogram.wcs.wcs_world2pix(world_coordinates, 0)

    lbv_pixels = iv.dendrogram.wcs.wcs_world2pix(world_coordinates, 0)

    l_lbv_pixels = lbv_pixels[:,0]
    b_lbv_pixels = lbv_pixels[:,1]
    v_lbv_pixels = lbv_pixels[:,2]

    l_scale_lbv = 0.125
    b_scale_lbv = 0.125
    v_scale_lbv = 1.3 

    lv_ells = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], 
                       width=2*garcia_sky_radius[i]/l_scale_lbv, 
                       height=garcia_catalog['DV'][i]/v_scale_lbv) for i in range(len(garcia_catalog))]
    for e in lv_ells:
        iv.ax_lv.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor('blue')
        e.set_zorder(2)

    iv.fig.canvas.draw()

