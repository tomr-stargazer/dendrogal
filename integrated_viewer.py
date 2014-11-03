"""
A viewer for astrodendro that looks like an integrated l,b and l,v thing.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from wcsaxes import WCSAxes

class IntegratedViewer(object):
    def __init__(self, dendrogram, hub, alignment='horizontal', cmap=plt.cm.gray,
                 clip_velocity=None, aspect=2.5):

        if dendrogram.data.ndim != 3:
            raise ValueError(
                "Only 3-dimensional arrays are supported")

        if not hasattr(hub, 'select_subtree'):
            raise NotImplementedError("astrodendro does not have scatter_picker enabled")

        self.hub = hub

        self.datacube = dendrogram.data
        self.dendrogram = dendrogram

        self.selected_contours = {} # selection_id -> (contour, contour) tuple

        if alignment == 'horizontal':
            figsize = (10, 4.4)
            ax_lb_limits = [0.1, 0.05, 0.8, 0.4]
            ax_lv_limits = [0.1, 0.5, 0.8, 0.5]
        elif alignment == 'vertical':
            figsize = (8, 6)
            ax_lb_limits = [0.1, 0.1, 0.375, 0.8]
            ax_lv_limits = [0.55, 0.1, 0.375, 0.8]
        elif alignment == 'pp' or alignment == 'lb':
            figsize = (10, 4)
            ax_lb_limits = [0.1, 0.1, 0.8, 0.8]
            ax_lv_limits = [0,0,0.01,0.01]
        else:
            raise ValueError("`alignment` must be 'horizontal' or 'vertical'")

        self.fig = plt.figure(figsize=figsize)
        self.cmap = cmap

        if self.dendrogram.wcs is not None:
            ax_lb = WCSAxes(self.fig, ax_lb_limits, wcs=self.dendrogram.wcs, slices=('x', 'y', 0))
            self.ax_lb = self.fig.add_axes(ax_lb)

            ax_lv = WCSAxes(self.fig, ax_lv_limits, wcs=self.dendrogram.wcs, slices=('x', 0, 'y'))
            self.ax_lv = self.fig.add_axes(ax_lv)
        else:
            self.ax_lb = self.fig.add_axes(ax_lb_limits)
            self.ax_lv = self.fig.add_axes(ax_lv_limits)

        array_lb = np.nansum(self.datacube, axis=0)
        array_lb[(array_lb < 0) | np.isinf(array_lb) | np.isnan(array_lb)] = 0
        self.array_lb = array_lb

        array_lv = np.nansum(self.datacube, axis=1)
        array_lv[(array_lv < 0) | np.isinf(array_lv) | np.isnan(array_lv)] = 0
        self.array_lv = array_lv

        if clip_velocity is None:
            if np.shape(array_lv)[0]*2.5 > np.shape(array_lb)[0]:
                self.clip_velocity = True
            else:
                self.clip_velocity = False
        else:
            self.clip_velocity = clip_velocity

        self.aspect = aspect

        self._draw_plot()
        self.hub.add_callback(self.update_selection)
        self.fig.canvas.mpl_connect('button_press_event', self.select_from_map)

        # If things are already selected in the hub, go select them!
        for selection_id in self.hub.selections:
            self.update_selection(selection_id)

    def _draw_plot(self):
    	""" Create an image to plot things onto. """

        len_v, len_l, len_b = self.datacube.shape

        self.image_lb = self.ax_lb.imshow(np.log10(self.array_lb+1), origin='lower', 
            interpolation='nearest', cmap=self.cmap)

        self.image_lv = self.ax_lv.imshow(np.log10(self.array_lv+1), origin='lower', 
            interpolation='nearest', cmap=self.cmap, aspect=self.aspect)

        # Trim the top and bottom of the l, v plot for cosmetic reasons
        if self.clip_velocity:
            crop_factor = np.round(len_v/5)
            self.ax_lv.set_ylim(crop_factor, len_v - crop_factor)

        self.fig.canvas.draw()

    def remove_contour(self, selection_id):
        if selection_id in self.selected_contours:
            for ax, contour in zip([self.ax_lb, self.ax_lv], self.selected_contours[selection_id]):
                for collection in contour.collections:
                    ax.collections.remove(collection)
            del self.selected_contours[selection_id]

    def select_from_map(self, event):

        # Only do this if no tools are currently selected
        if event.canvas.toolbar.mode != '':
            return
        if event.button not in self.hub.colors:
            return

        input_key = event.button

        if event.inaxes is self.ax_lb:

            # Find pixel co-ordinates of click
            ix = int(round(event.xdata))
            iy = int(round(event.ydata))

            try:
                iz = np.where(self.datacube[:, iy, ix] == np.nanmax(self.datacube[:, iy, ix]))[0][0]
            except IndexError:
                iz = 0

        elif event.inaxes is self.ax_lv:

            # Find pixel co-ordinates of click
            ix = int(round(event.xdata))
            iz = int(round(event.ydata))

            try:
                iy = np.where(self.datacube[iz, :, ix] == np.nanmax(self.datacube[iz, :, ix]))[0][0]
            except IndexError:
                iy = 0

        else:
            return

        indices = (iz, iy, ix)

        # Select the structure
        structure = self.dendrogram.structure_at(indices)
        self.hub.select(input_key, structure)
        self.update_selection(input_key)

        # Re-draw
        event.canvas.draw()

    def update_selection(self, selection_id):
        """Highlight seleted structures"""

        self.remove_contour(selection_id)
        
        if selection_id in self.hub.selections:

            structures = self.hub.selections[selection_id]
            select_subtree = self.hub.select_subtree[selection_id]

            struct = structures[0]

            if struct is not None:
                if select_subtree:
                    mask = struct.get_mask(subtree=True)
                else:
                    mask = reduce(np.add, [structure.get_mask(subtree=True) for structure in structures])
                # numpy magic sets nonzero values to True
                mask_lb = np.sum(mask, axis=0).astype('bool') 
                mask_lv = np.sum(mask, axis=1).astype('bool')

                self.selected_contours[selection_id] = (
                    self.ax_lb.contour(
                        mask_lb, colors=self.hub.colors[selection_id],
                        linewidths=5, levels=[0.5], alpha=0.9, zorder=struct.height), 
                    self.ax_lv.contour(
                        mask_lv, colors=self.hub.colors[selection_id],
                        linewidths=5, levels=[0.5], alpha=0.9, zorder=struct.height) )

        self.fig.canvas.draw()
