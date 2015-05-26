"""
Makes a thumbnail showing a cloud in data-space and dendro-space.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from dendrogal.production.integrated_map_figure import integrated_map_axes_lb, integrated_map_axes_lv

path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'


def single_cloud_lb_thumbnail():
    """
    Makes an l, b thumbnail of a single cloud.

    Shows a contour of the emission region used, and a representative ellipse.

    Integrates only over the relevant velocities (twice the RMS in each direction).

    """
    pass


def single_cloud_lv_thumbnail():
    """
    Makes an l, v thumbnail of a single cloud.

    Shows a contour of the emission region used, and a representative ellipse.

    Integrates only over the relevant latitudes (twice the sigma_r in each direction).

    """    
    pass


def single_cloud_dendro_thumbnail(ax, dendrogram, cloud_idx):
    """
    Makes a dendrogram thumbnail for a single cloud.

    Highlights the appropriate branches (and substructure of course).

    Zooms the dendrogram to the right location.

    """

    d = dendrogram

    p = d.plotter()

    p.plot_tree(ax, color='black')

    p.plot_tree(ax, structure=[d[cloud_idx]], color=colorbrewer_blue, lw=2)

    # everything below is machinery to zoom the figure reasonably
    structures = d[cloud_idx].descendants + [d[cloud_idx]]
    structure_positions = [p._cached_positions[x] for x in structures]
    min_position = min(structure_positions)
    max_position = max(structure_positions)
    range_positions = max_position - min_position

    ax.set_xlim(min_position - range_positions/4, max_position + range_positions/4)

    vmin_list = [x.vmin for x in structures]
    vmax_list = [x.vmax for x in structures]

    min_vmin = min(vmin_list)
    max_vmax = max(vmax_list)
    v_range = max_vmax - min_vmin

    ax.set_ylim(min_vmin - v_range/10, max_vmax + v_range/10)

    return p



def make_thumbnail_dendro_figure(dendrogram, catalog, cloud_idx):

    d = dendrogram
    data = d.data
    struct = d[cloud_idx]

    cloud_row = catalog[catalog['_idx'] == cloud_idx]

    velocity_integration = [cloud_row['v_cen']-2*cloud_row['v_rms'],
                            cloud_row['v_cen']+2*cloud_row['v_rms']]

    print velocity_integration

    fig = plt.figure()

    ax_dendro = fig.add_subplot(122)

    single_cloud_dendro_thumbnail(ax_dendro, d, cloud_idx)

    # draw a map on ax_map
    ax_map_limits = [0.1, 0.55, 0.35, 0.35]
    ax_map = integrated_map_axes_lb(fig, ax_map_limits, data, d.wcs)

    ax_lv_limits =  [0.1, 0.1, 0.35, 0.35]  
    ax_lv = integrated_map_axes_lv(fig, ax_lv_limits, data, d.wcs)      

    # A. draw ellipse
    world_coordinates = np.vstack([cloud_row['x_cen'], cloud_row['y_cen'], cloud_row['v_cen']]).T
    lbv_pixels = d.wcs.wcs_world2pix(world_coordinates, 0)

    l_lbv_pixels = lbv_pixels[:,0]
    b_lbv_pixels = lbv_pixels[:,1]
    v_lbv_pixels = lbv_pixels[:,2]

    l_scale_lbv = ax_lv.spatial_scale.value
    b_scale_lbv = ax_lv.spatial_scale.value
    v_scale_lbv = ax_lv.velocity_scale.value

    lb_ell = [Ellipse(xy=zip(l_lbv_pixels, b_lbv_pixels)[i], 
                          angle=cloud_row['position_angle'][i],
                          width=2*cloud_row['major_sigma'][i]/l_scale_lbv, 
                          height=2*cloud_row['minor_sigma'][i]/b_scale_lbv) 
                  for i in range(len(cloud_row))]

    lb_ell2 = [Ellipse(xy=zip(l_lbv_pixels, b_lbv_pixels)[i], 
                          angle=cloud_row['position_angle'][i],
                          width=2*cloud_row['major_sigma'][i]/l_scale_lbv, 
                          height=2*cloud_row['minor_sigma'][i]/b_scale_lbv) 
                  for i in range(len(cloud_row))]

    for e in lb_ell:
        ax_map.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor(colorbrewer_red)
        e.set_linewidth(1.25)
        e.set_zorder(0.95)

    for e in lb_ell2:
        ax_map.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor('white')
        e.set_linewidth(4)
        e.set_alpha(0.8)
        e.set_zorder(0.9)

    lv_ells = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], 
                       width=2*cloud_row['major_sigma'][i]/l_scale_lbv, 
                       height=2*cloud_row['v_rms'][i]/v_scale_lbv) for i in range(len(cloud_row))]

    lv_ells2 = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], 
                       width=2*cloud_row['major_sigma'][i]/l_scale_lbv, 
                       height=2*cloud_row['v_rms'][i]/v_scale_lbv) for i in range(len(cloud_row))]

    for e in lv_ells:
        ax_lv.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor(colorbrewer_red)
        e.set_linewidth(1.25)
        e.set_zorder(0.95)

    for e in lv_ells2:
        ax_lv.add_artist(e)
        e.set_facecolor('none')
        e.set_edgecolor('white')
        e.set_linewidth(4)
        e.set_alpha(0.8)
        e.set_zorder(0.9)


    # draw a dendrogram on ax_dendro

