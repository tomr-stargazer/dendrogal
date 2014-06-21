"""
This is a test to see if I can overlay some galactic coordinates 
and scatter plotted points on a Robert Hurt galaxy image.

TECHNICALLY most of these are actually heliocentric XYZ coordinates 
oriented w.r.t. the Galaxy's center and plane.

"""

from __future__ import division
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imread

from astropy import units as u
from astropy import constants as c

dropbox_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/")

hurt_image_location = dropbox_path+"hurt_rotated_small.jpg.pagespeed.ce.PXjylj-jgc.jpg"
hurt_image = imread(hurt_image_location)

# These numbers were derived using WebPlotDigitizer
# assuming that the galactic center is 8340 pc from the Sun
# (i.e., disregarding the distance axes on the image itself!)

px_between_sun_and_galactic_center = 115.8
distance_between_sun_and_galactic_center = 8340 # pc

x_len, y_len, a = np.shape(hurt_image)

image_length_in_pc = x_len * distance_between_sun_and_galactic_center/px_between_sun_and_galactic_center

galactic_x_extent = [image_length_in_pc/2 + distance_between_sun_and_galactic_center,
	                 -image_length_in_pc/2 + distance_between_sun_and_galactic_center]
galactic_y_extent = [-image_length_in_pc/2, image_length_in_pc/2]

# galactic_y_extent = [-18378.6, 18269.7] # "y" axis runs through 90 degrees l and 270 degrees l
# galactic_x_extent = [25057.8, -11329.5] # "x" axis runs through Galactic Center and Anticenter

fig = plt.figure()
plt.imshow(hurt_image, extent=galactic_y_extent+galactic_x_extent)

plt.scatter(np.arange(100), 100*np.sin(np.arange(100)/10), marker='o', c='red')
plt.show()