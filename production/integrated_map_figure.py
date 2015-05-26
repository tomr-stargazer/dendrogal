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

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'

def latitude_world2pix(wcs, latitudes):
    # assumes an l, b, v datacube (or at least that axis=1 is the latitude axis)

    filler = np.repeat(0, len(latitudes))

    world_coordinates = np.vstack([filler, latitudes, filler]).T
    pixels = wcs.wcs_world2pix(world_coordinates, 0)

    return pixels[:,1]


def velocity_world2pix(wcs, velocities):
    # assumes an l, b, v datacube (or at least that axis=2 is the velocity axis)

    filler = np.repeat(0, len(velocities))

    world_coordinates = np.vstack([filler, filler, velocities]).T
    pixels = wcs.wcs_world2pix(world_coordinates, 0)

    return pixels[:,2]


class IntegratedMapFigure(object):

    def __init__(self, datacube, wcs_object,  
                 aspect_in_units=5*(u.km/u.s)/u.deg, cmap='RdGy_r',
                 integration_limits=(-1,1), figsize=(10,9) ):

        self.datacube = datacube
        self.wcs = wcs_object

        self.fig = plt.figure(figsize=figsize)
        self.cmap = cmap

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

        self.fig.canvas.draw()


    def draw_contours(self, mask_3d, **kwargs):
        """
        Draws contours when you pass in a mask array!
        """

        mask_lv = np.sum(mask_3d, axis=1).astype('bool')

        self.ax.contour(mask_lv, **kwargs)

        self.fig.canvas.draw()


def integrated_map_axes_lb(fig, ax_limits, datacube, wcs_object, integration_limits,
                           cmap='RdGy_r' ):
    """ 
    Hackish way of exporting the above behavior to an ax object.

    But it's actually on-the-sky.

    """

    ax = WCSAxes(fig, ax_limits, wcs=wcs_object, 
                     slices=('x', 'y', 0))
    fig.add_axes(ax)

    min_v_px, max_v_px = velocity_world2pix(wcs_object, integration_limits)

    array_lb = np.nansum(datacube[min_v_px:max_v_px,:,:], axis=0)
    array_lb[(array_lb < 0) | np.isinf(array_lb) | np.isnan(array_lb)] = 0

    ax.spatial_scale = np.abs(wcs_object.wcs.cdelt[0]) * u.deg
    ax.velocity_scale = np.abs(wcs_object.wcs.cdelt[2]) * u.km/u.s    

    image_lb = ax.imshow(np.log10(array_lb+1), origin='lower', 
        interpolation='nearest', cmap=cmap, aspect=1)

    return ax


def integrated_map_axes_lv(fig, ax_limits, datacube, wcs_object, integration_limits,
                           aspect_in_units=5*(u.km/u.s)/u.deg,
                           cmap='RdGy_r'):
    """ 
    Hackish way of exporting the above behavior to an ax object.

    """

    ax = WCSAxes(fig, ax_limits, wcs=wcs_object, 
                     slices=('x', 0, 'y'))
    fig.add_axes(ax)

    min_b_px, max_b_px = latitude_world2pix(wcs_object, integration_limits)

    array_lv = np.nansum(datacube[:,min_b_px:max_b_px,:], axis=1)
    array_lv[(array_lv < 0) | np.isinf(array_lv) | np.isnan(array_lv)] = 0

    ax.spatial_scale = np.abs(wcs_object.wcs.cdelt[0]) * u.deg
    ax.velocity_scale = np.abs(wcs_object.wcs.cdelt[2]) * u.km/u.s

    ax.aspect_px = (ax.velocity_scale/ax.spatial_scale)/aspect_in_units

    image_lv = ax.imshow(np.log10(array_lv+1), origin='lower', 
        interpolation='nearest', cmap=cmap, aspect=ax.aspect_px)

    return ax    

