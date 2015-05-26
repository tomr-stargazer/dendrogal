"""
Makes a thumbnail showing a cloud in data-space and dendro-space.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from dendrogal.production.integrated_map_figure import integrated_map_axes_lb

path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'


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

    # draw a map on ax_map
    ax_map_limits = [0.1, 0.1, 0.35, 0.8]
    ax_map = integrated_map_axes_lb(fig, ax_map_limits, data, d.wcs)

    # A. draw ellipse
    world_coordinates = np.vstack([cloud_row['x_cen'], cloud_row['y_cen'], cloud_row['v_cen']]).T
    lbv_pixels = d.wcs.wcs_world2pix(world_coordinates, 0)

    l_lbv_pixels = lbv_pixels[:,0]
    b_lbv_pixels = lbv_pixels[:,1]

    l_scale_lbv = 0.125
    b_scale_lbv = 0.125

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


    # draw a dendrogram on ax_dendro

