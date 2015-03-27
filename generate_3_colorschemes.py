"""
Script that makes 3 colorschemes for the integrated map / ellipse figures.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from astrodendro_analysis.production.map_figures import make_quadrant_lbv_map
from dummy_generator_of_d_and_catalog import d, cloud_catalog

import alt_cmapper
dissapgr_cmap = alt_cmapper.makecmap("dissapgr", rev=False, ncol=50)
alarm_cmap = alt_cmapper.makecmap("alarm.p2.0.5", rev=False, ncol=50)
# here's the good stuff

cmap_list = ['RdGy_r', dissapgr_cmap, alarm_cmap]

for cmap in cmap_list:

    make_quadrant_lbv_map(cloud_catalog, d, contour_select=False, aspect=1, cmap=cmap)

# def make_quadrant_lbv_map(cloud_catalog, dendgrogram, contour_select=True, 
#                           alignment='lv', vscale=0.65, xscale=0.125, 
#                           aspect=1/3, cmap='gray_r', **kwargs):

plt.show()