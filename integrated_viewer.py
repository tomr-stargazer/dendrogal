"""
A viewer for astrodendro that looks like an integrated l,b and l,v thing.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from wcsaxes import WCSAxes

class IntegratedViewer(object):
    def __init__(self, dendrogram, hub, wcs=None):

        if dendrogram.data.ndim != 3:
            raise ValueError(
                "Only 3-dimensional arrays are supported")

        if not hasattr(hub, 'select_subtree'):
            raise NotImplementedError("astrodendro does not have scatter_picker enabled")

        self.hub = hub
        self.hub.add_callback(self.update_selection)

        self.datacube = dendrogram.data
        self.dendrogram = dendrogram

        self.selected_contours = {} # selection_id -> (contour, contour) tuple

        self.fig = plt.figure(figsize=(10, 4.4))

        ax_lb_limits = [0.1, 0.05, 0.8, 0.4]
        ax_lv_limits = [0.1, 0.5, 0.8, 0.5]

        if wcs is not None:
            ax_lb = WCSAxes(self.fig, ax_lb_limits, wcs=wcs)
            self.ax_lb = self.fig.add_axes(ax_lb)
            # ax_lv = WCSAxes(self.fig, self.ax_lv_limits, wcs=wcs)
            # self.ax_lv = self.fig.add_axes(ax_lv)
            self.ax_lv = self.fig.add_axes(ax_lv_limits)             #temporary
        else:
            self.ax_lb = self.fig.add_axes(ax_lb_limits)
            self.ax_lv = self.fig.add_axes(ax_lv_limits)

        self.array_lb = np.nansum(self.datacube, axis=0)
        self.array_lv = np.nansum(self.datacube, axis=1)

        self._draw_plot()
        self.hub.add_callback(self.update_selection)
        self.fig.canvas.mpl_connect('button_press_event', self.select_from_map)


    def _draw_plot(self):
    	""" Create an image to plot things onto. """

        len_v, len_l, len_b = self.datacube.shape

        self._clim_lb = (np.min(self.array_lb[~np.isnan(self.array_lb) & ~np.isinf(self.array_lb)]),
        	             np.max(self.array_lb[~np.isnan(self.array_lb) & ~np.isinf(self.array_lb)]))

        self._clim_lv = (np.min(self.array_lv[~np.isnan(self.array_lv) & ~np.isinf(self.array_lv)]),
        	             np.max(self.array_lv[~np.isnan(self.array_lv) & ~np.isinf(self.array_lv)]))

        self.image_lb = self.ax_lb.imshow(self.array_lb, origin='lower', 
            interpolation='nearest', vmin=self._clim_lb[0], 
            vmax=0.025*self._clim_lb[1], cmap=plt.cm.gray)

        self.image_lv = self.ax_lv.imshow(self.array_lv, origin='lower', 
            interpolation='nearest', vmin=self._clim_lv[0], 
            vmax=0.05*self._clim_lv[1], cmap=plt.cm.gray, aspect=2.5)

        # Trim the top and bottom of the l, v plot for cosmetic reasons
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

            iz = np.where(self.datacube[:, iy, ix] == np.nanmax(self.datacube[:, iy, ix]))[0][0]

        elif event.inaxes is self.ax_lv:

            # Find pixel co-ordinates of click
            ix = int(round(event.xdata))
            iz = int(round(event.ydata))

            iy = np.where(self.datacube[iz, :, ix] == np.nanmax(self.datacube[iz, :, ix]))[0][0]

        else:
            return

        indices = (iz, iy, ix)
        assert type(iz) == type(iy) == type(ix) == type(1)

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
                        linewidths=2, levels=[0.5], alpha=0.75, zorder=struct.height), 
                    self.ax_lv.contour(
                        mask_lv, colors=self.hub.colors[selection_id],
                        linewidths=2, levels=[0.5], alpha=0.75, zorder=struct.height) )

        self.fig.canvas.draw()
