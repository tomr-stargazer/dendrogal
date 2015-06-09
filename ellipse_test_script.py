""" A cute script to test which ellipse parameters work best, visually, on the RdGy_r cmap. """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from dummy_generator_of_d_and_catalog import d, cloud_catalog

from dendrogal.production.integrated_map_figure import IntegratedMapFigure

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'

imf = IntegratedMapFigure(d.data, d.wcs, figsize=(10,9))

v_centers = np.linspace(-75, 100, 10) + 10
l_centers = np.ones_like(v_centers) * 35
b_centers = np.zeros_like(v_centers)

world_coordinates = np.vstack([l_centers, b_centers, v_centers]).T

lbv_pixels = d.wcs.wcs_world2pix(world_coordinates, 0)

l_lbv_pixels = lbv_pixels[:,0]
b_lbv_pixels = lbv_pixels[:,1]
v_lbv_pixels = lbv_pixels[:,2]

lv_ells = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], width=10, height=10) for i in range(len(v_centers))]
lv_ells2 = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], width=10, height=10) for i in range(len(v_centers))]
lv_ells3 = [Ellipse(xy=zip(l_lbv_pixels, v_lbv_pixels)[i], width=10, height=10) for i in range(len(v_centers))]

# for e in lv_ells:
#     imf.ax.add_artist(e)
#     e.set_facecolor('none')
#     e.set_edgecolor('black')
#     e.set_linewidth(1)
#     e.set_zorder(1)

for e in lv_ells:
    imf.ax.add_artist(e)
    e.set_facecolor('none')
    e.set_edgecolor(colorbrewer_red)
    e.set_linewidth(2)
    e.set_zorder(0.95)

for e in lv_ells2:
    imf.ax.add_artist(e)
    e.set_facecolor('none')
    e.set_edgecolor('white')
    e.set_linewidth(4)
    e.set_alpha(0.8)
    e.set_zorder(0.9)

imf.fig.canvas.draw()
plt.show()