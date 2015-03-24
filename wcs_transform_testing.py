""" for testing wcsaxes transforms """

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from numpy.testing import assert_allclose, assert_equal

from astropy import wcs
import wcsaxes

header = (""
    "WCSAXES =                    3 / Number of coordinate axes                      "
    "CRPIX1  =                  1.0 / Pixel coordinate of reference point            "
    "CRPIX2  =                 33.0 / Pixel coordinate of reference point            "
    "CRPIX3  =                  1.0 / Pixel coordinate of reference point            "
    "CDELT1  =                -0.25 / [deg] Coordinate increment at reference point  "
    "CDELT2  =                 0.25 / [deg] Coordinate increment at reference point  "
    "CDELT3  =              1.30038 / [m/s] Coordinate increment at reference point  "
    "CUNIT1  = 'deg'                / Units of coordinate increment and value        "
    "CUNIT2  = 'deg'                / Units of coordinate increment and value        "
    "CUNIT3  = 'm/s'                / Units of coordinate increment and value        "
    "CTYPE1  = 'GLON-CAR'           / galactic longitude, plate caree projection     "
    "CTYPE2  = 'GLAT-CAR'           / galactic latitude, plate caree projection      "
    "CTYPE3  = 'VOPT'               / Optical velocity (linear)                      "
    "CRVAL1  =                 75.0 / [deg] Coordinate value at reference point      "
    "CRVAL2  =                  0.0 / [deg] Coordinate value at reference point      "
    "CRVAL3  =             -90.2712 / [m/s] Coordinate value at reference point      "
    "LONPOLE =                  0.0 / [deg] Native longitude of celestial pole       "
    "LATPOLE =                 90.0 / [deg] Native latitude of celestial pole        "
    "SPECSYS = 'LSRK'               / Reference frame of spectral coordinates        END")

wcs_object = wcs.wcs.WCS(header)

# >>> wcs_object.printwcs()
# WCS Keywords

# Number of WCS axes: 3
# CTYPE : 'GLON-CAR'  'GLAT-CAR'  'VOPT'
# CRVAL : 75.0  0.0  -90.271199999999993
# CRPIX : 1.0  33.0  1.0
# PC1_1 PC1_2 PC1_3 : 1.0  0.0  0.0
# PC2_1 PC2_2 PC2_3 : 0.0  1.0  0.0
# PC3_1 PC3_2 PC3_3 : 0.0  0.0  1.0
# CDELT : -0.25  0.25  1.3003800000000001
# NAXIS    : 0 0

world_coordinates = np.vstack([np.arange(10), np.repeat(0, 10), np.arange(10)]).T

pixel_from_wcs = wcs_object.wcs_world2pix(world_coordinates, 0)

x_from_wcs = pixel_from_wcs[:,0]
y_from_wcs = pixel_from_wcs[:,2]

fig = plt.figure()

ax = wcsaxes.WCSAxes(fig, [0.1 ,0.1 ,0.8, 0.8], wcs=wcs_object, slices=('x', 0, 'y'))
fig.add_axes(ax)

# ax.plot(np.arange(10), np.arange(10), transform=ax.get_transform('world'))
tr = ax.get_transform('world')
pixel_from_tr = tr.transform(np.array([np.arange(10),
                               np.repeat(-90, 10),
                               np.arange(10)]).transpose())
x_from_tr, y_from_tr = pixel_from_tr[:,0], pixel_from_tr[:,1]
ax.plot(x_from_tr, y_from_tr, 'b')
ax.plot(x_from_wcs, y_from_wcs, 'r')
plt.show()
fig.canvas.draw()

assert_equal(x_from_tr, x_from_wcs)
assert_equal(y_from_tr, y_from_wcs)
