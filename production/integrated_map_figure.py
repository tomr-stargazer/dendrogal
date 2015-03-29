"""
A class that makes Integrated Map Figures. 

Inspired largely by IntegratedViewer, but 
(a) with the interactive/synchronization with a Viewer.hub bits stripped out
(b) more fine control over integration ranges
(c) better use of WCS information.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from wcsaxes import WCSAxes

def latitude_world2pix(wcs, latitudes):
    # assumes an l, b, v datacube (or at least that axis=1 is the latitude axis)

    filler = np.repeat(0, len(latitudes))

    world_coordinates = np.vstack([filler, latitudes, filler]).T
    pixels = wcs.wcs_world2pix(world_coordinates, 0)

    return pixels[:,1]

class IntegratedMapFigure(object):

    def __init__(self, datacube, wcs_object, contour=None, 
                 aspect_in_units=5*(u.km/u.s)/u.deg, cmap='RdGy_r',
                 integration_limits=(-1,1), figsize=(10,6) ):

        self.datacube = datacube
        self.wcs = wcs_object

        self.fig = plt.figure(figsize=figsize)
        self.cmap = cmap

        self.contour = contour

        ax = WCSAxes(self.fig, [0.1, 0.1, 0.8, 0.8], wcs=self.wcs, 
                     slices=('x', 0, 'y'))
        self.ax = self.fig.add_axes(ax)

        min_b_px, max_b_px = latitude_world2pix(self.wcs, integration_limits)

        array_lv = np.nansum(self.datacube[:,min_b_px:max_b_px,:], axis=1)
        array_lv[(array_lv < 0) | np.isinf(array_lv) | np.isnan(array_lv)] = 0
        self.array_lv = array_lv

        self.spatial_scale = np.abs(self.wcs.wcs.cdelt[0]) * u.deg
        self.velocity_scale = np.abs(self.wcs.wcs.cdelt[2]) * u.km/u.s

        self.aspect_px = (self.velocity_scale/self.spatial_scale)/aspect_in_units

        self._draw_plot()

    def _draw_plot(self):

        self.image_lv = self.ax.imshow(np.log10(self.array_lv+1), origin='lower', 
            interpolation='nearest', cmap=self.cmap, aspect=self.aspect_px)


