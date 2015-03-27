"""
Script that makes 3 colorschemes for the integrated map / ellipse figures.

"""

from __future__ import division

import os

import numpy as np
import matplotlib.pyplot as plt

from astrodendro_analysis.production.map_figures import make_quadrant_lbv_map
from dummy_generator_of_d_and_catalog import d, cloud_catalog

save_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/color_schemes/")

import alt_cmapper
butterflyfairy_cmap = alt_cmapper.makecmap("butterflyfairy", rev=True, ncol=50)
alarm_cmap = alt_cmapper.makecmap("alarm.p2.0.5", rev=False, ncol=50)
# here's the good stuff

cmap_list = ['gray_r', 'spectral', 'RdGy_r', butterflyfairy_cmap, alarm_cmap]
cmap_name_list = ['gray_r', 'spectral', 'RdGy_r', 'butterflyfairy', 'alarm']

for cmap, name in zip(cmap_list, cmap_name_list):

    iv = make_quadrant_lbv_map(cloud_catalog, d, contour_select=False, aspect=1, cmap=cmap)

    iv.fig.set_figwidth(8.7)
    iv.fig.set_figheight(7.3)

    iv.fig.savefig(save_path+"{0}_map.png".format(name))

plt.show()