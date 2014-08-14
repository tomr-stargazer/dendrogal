"""
This is a script that helps me overlay the Dame (2001) "cartoon" over integrated maps.

"""

from __future__ import division
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imread

from wcsaxes import WCSAxes

dropbox_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/")

dame_lb_cartoon_location = dropbox_path+"Fig2_Dame_cropped.png"
lb_cartoon_image = imread(dame_lb_cartoon_location)

test_fig = plt.figure()
test_ax = test_fig.add_subplot(111)

def overlay_lb_cartoon(ax, wcs):


	

	# 	plt.imshow(hurt_image, zorder=-10, extent=extents, **kwargs)
	

	pass


