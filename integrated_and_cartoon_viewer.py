"""
A little moduel that creates the IntegratedViewer knockoff that has a Dame cartoon.

The axes of the two views are shared.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from wcsaxes import WCSAxes

from overlay_dame_cartoon import overlay_lb_cartoon

class CartoonDualViewer(object):
    def __init__(self, dendrogram, hub, alignment='horizontal', cmap=plt.cm.gray,
                 clip_velocity=None):

        if dendrogram.data.ndim != 3:
            raise ValueError(
                "Only 3-dimensional arrays are supported")

        if not hasattr(hub, 'select_subtree'):
            raise NotImplementedError("astrodendro does not have scatter_picker enabled")

        self.hub = hub

        self.datacube = dendrogram.data
        self.dendrogram = dendrogram
        self.wcs = dendrogram.wcs

        self.selected_contours = {} # selection_id -> (contour, contour) tuple

        if alignment == 'horizontal':
            figsize = (11, 9)
            ax_integrated_limits = [0.05, 0.05, 0.9, 0.5]
            ax_cartoon_limits = [0.05, 0.55, 0.9, 0.5]
        elif alignment == 'vertical':
            figsize = (8, 6)
            ax_integrated_limits = [0.1, 0.1, 0.375, 0.8]
            ax_cartoon_limits = [0.55, 0.1, 0.375, 0.8]
        # elif alignment == 'pp' or alignment == 'lb':
        #     figsize = (10, 4)
        #     ax_integrated_limits = [0.1, 0.1, 0.8, 0.8]
        #     ax_cartoon_limits = [0,0,0.01,0.01]
        else:
            raise ValueError("`alignment` must be 'horizontal' or 'vertical'")

        self.fig = plt.figure(figsize=figsize)
        self.cmap = cmap

        if self.dendrogram.wcs is not None:
            ax_integrated = WCSAxes(self.fig, ax_integrated_limits, wcs=self.dendrogram.wcs, slices=('x', 'y', 0))
            self.ax_integrated = self.fig.add_axes(ax_integrated)

            ax_cartoon = WCSAxes(self.fig, ax_cartoon_limits, wcs=self.dendrogram.wcs, slices=('x', 'y', 0), 
            	sharex=ax_integrated, sharey=ax_integrated)
            self.ax_cartoon = self.fig.add_axes(ax_cartoon)
        else:
            self.ax_integrated = self.fig.add_axes(ax_integrated_limits)
            self.ax_cartoon = self.fig.add_axes(ax_cartoon_limits,
            	sharex=ax_integrated, sharey=ax_integrated)

        array_lb = np.nansum(self.datacube, axis=0)
        array_lb[(array_lb < 0) | np.isinf(array_lb) | np.isnan(array_lb)] = 0
        self.array_lb = array_lb

        # array_lv = np.nansum(self.datacube, axis=1)
        # array_lv[(array_lv < 0) | np.isinf(array_lv) | np.isnan(array_lv)] = 0
        # self.array_lv = array_lv

        # if clip_velocity is None:
        #     if np.shape(array_lv)[0]*2.5 > np.shape(array_lb)[0]:
        #         self.clip_velocity = True
        #     else:
        #         self.clip_velocity = False
        # else:
        #     self.clip_velocity = clip_velocity

        self._draw_plot()
        self.hub.add_callback(self.update_selection)
        self.fig.canvas.mpl_connect('button_press_event', self.select_from_map)

        # If things are already selected in the hub, go select them!
        for selection_id in self.hub.selections:
            self.update_selection(selection_id)

    def _draw_plot(self):
    	""" Create an image to plot things onto. """

        len_v, len_l, len_b = self.datacube.shape

        self.image_integrated = self.ax_integrated.imshow(np.log10(self.array_lb+1), origin='lower', 
            interpolation='nearest', cmap=self.cmap)

        self.image_cartoon = overlay_lb_cartoon(self.ax_cartoon)[0]

        # Trim the top and bottom of the l, v plot for cosmetic reasons
        # if self.clip_velocity:
        #     crop_factor = np.round(len_v/5)
        #     self.ax_cartoon.set_ylim(crop_factor, len_v - crop_factor)

        self.fig.canvas.draw()

    def remove_contour(self, selection_id):
        if selection_id in self.selected_contours:
            for ax, contour in zip([self.ax_integrated, self.ax_cartoon], self.selected_contours[selection_id]):
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

        if event.inaxes is self.ax_integrated or event.inaxes is self.ax_cartoon:

            # Find pixel co-ordinates of click
            ix = int(round(event.xdata))
            iy = int(round(event.ydata))

            try:
                iz = np.where(self.datacube[:, iy, ix] == np.nanmax(self.datacube[:, iy, ix]))[0][0]
            except IndexError:
                iz = 0
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
                    self.ax_integrated.contour(
                        mask_lb, colors=self.hub.colors[selection_id],
                        linewidths=5, levels=[0.5], alpha=0.9, zorder=struct.height), 
                    self.ax_cartoon.contour(
                        mask_lb, colors=self.hub.colors[selection_id],
                        linewidths=5, levels=[0.5], alpha=0.9, zorder=struct.height) )

        self.fig.canvas.draw()
