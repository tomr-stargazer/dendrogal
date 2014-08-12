"""
Snippet of plotting code that reproduces my WCSAxes glitch.

"""

from __future__ import division

import os.path
import matplotlib.pyplot as plt

from astropy import wcs
from astropy.io.fits import getdata

from wcsaxes import WCSAxes
from demo import downsample_and_transpose_data_and_header

data_path = os.path.expanduser("~/Dropbox/College/Astro99/DATA/")

data_file = 'COGAL_local_mom.fits'

datacube, datacube_header = getdata(data_path+data_file, memmap=True,
	header=True)

datacube_dt, datacube_dt_header = \
	downsample_and_transpose_data_and_header(datacube, datacube_header, 1, (2,0,1), resample=False, recenter=True)
datacube_dt_wcs = wcs.wcs.WCS(datacube_dt_header)

datacube_dt_wcs.wcs.bounds_check(pix2world=False)

fig = plt.figure(figsize=(13,3.5))
ax_image_limits = [0.1, 0.1, 0.8, 0.8]
ax_image = fig.add_axes(WCSAxes(fig, ax_image_limits, wcs=datacube_dt_wcs, slices=('x','y', 0)))

image = ax_image.imshow(datacube_dt[50, :,:], origin='lower', interpolation='nearest', cmap=plt.cm.gray)

fig.show()