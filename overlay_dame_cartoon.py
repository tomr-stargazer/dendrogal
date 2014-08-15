"""
This is a script that helps me overlay the Dame (2001) "cartoon" over integrated maps.

"""

from __future__ import division
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imread

from wcsaxes import WCSAxes
from astropy import wcs

dropbox_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/")

dame_lb_cartoon_location = dropbox_path+"Fig2_Dame_cropped.png"
lb_cartoon_image = imread(dame_lb_cartoon_location)

# test_fig = plt.figure()
# test_ax = test_fig.add_subplot(111)

# in coordinate values, where are the corners?
corners = [
	(206.50549450549448, -49.60396039603961),
	(206.50549450549448, 44.05940594059403),
	(-190.87912087912053, 44.05940594059389),
	(-190.87912087912008, -49.60396039603967)
	]

px_per_degree = 4

def overlay_lb_cartoon(ax, wcs_object=None, **kwargs):

	if wcs_object is None:
		try:
			wcs_object= ax.wcs.sub([wcs.WCSSUB_CELESTIAL])
			wcs_object.wcs.bounds_check(False,False)
		except AttributeError:
			print "`ax` has no `wcs` -- please provide one"
			return

	left, a = wcs_object.all_world2pix(206.50549450549448-180, 0, 0) 
	right, a = wcs_object.all_world2pix(-190.87912087912053+360, 0, 0)
	b, top = wcs_object.all_world2pix(0, 44.05940594059403, 0)
	b, bottom = wcs_object.all_world2pix(0, -49.60396039603961, 0)

	# experimental...
	left2, a = wcs_object.all_world2pix(0, 0, 0) 
	right2, a = wcs_object.all_world2pix(0, 0, 0)

	print left2,

	left2 -= 206.50549450549448*px_per_degree
	right2 += (-190.87912087912053+360)*px_per_degree

	print left2, right2

	# extents = [left, right, bottom, top]
	extents = [-106, 1483.6, bottom, top]
	
	ax.imshow(lb_cartoon_image, zorder=0.05, extent=extents, **kwargs)
	
	return extents, wcs_object	

