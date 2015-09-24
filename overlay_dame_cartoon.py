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

code_directory = os.path.dirname(__file__)
cartoon_path = os.path.join(code_directory, "images")

dame_lb_cartoon_location = os.path.join(cartoon_path,"Fig2_Dame_cropped.png")
lb_cartoon_image = imread(dame_lb_cartoon_location)

dame_lv_cartoon_location = os.path.join(cartoon_path,"Fig3_Dame_cropped.png")
lv_cartoon_image = imread(dame_lv_cartoon_location)

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

	# print left2,

	left2 -= 206.50549450549448*px_per_degree
	right2 += (-190.87912087912053+360)*px_per_degree

	# print left2, right2

	# extents = [left, right, bottom, top]
	extents = [-106, 1483.6, bottom, top]
	
	image = ax.imshow(lb_cartoon_image, zorder=0.05, extent=extents, **kwargs)
	
	return image, extents, wcs_object

def overlay_lv_cartoon(ax, wcs_object=None, **kwargs):

	# if wcs_object is None:
	# 	try:
	# 		wcs_object= ax.wcs.sub([wcs.WCSSUB_CELESTIAL])
	# 		wcs_object.wcs.bounds_check(False,False)
	# 	except AttributeError:
	# 		print "`ax` has no `wcs` -- please provide one"
	# 		return

	# left, a = wcs_object.all_world2pix(206.50549450549448-180, 0, 0) 
	# right, a = wcs_object.all_world2pix(-190.87912087912053+360, 0, 0)
	# b, top = wcs_object.all_world2pix(0, 44.05940594059403, 0)
	# b, bottom = wcs_object.all_world2pix(0, -49.60396039603961, 0)

	# experimental...
	# left2, a = wcs_object.all_world2pix(0, 0, 0) 
	# right2, a = wcs_object.all_world2pix(0, 0, 0)

	# print left2,

	# left2 -= 206.50549450549448*px_per_degree
	# right2 += (-190.87912087912053+360)*px_per_degree

	# print left2, right2

	# extents = [left, right, bottom, top]
	extents = [-112.92, 1479.88, 39.46, 199.92]
	
	image = ax.imshow(lv_cartoon_image, zorder=0.05, extent=extents, aspect=2.5, **kwargs)
	
	return image, extents, wcs_object

