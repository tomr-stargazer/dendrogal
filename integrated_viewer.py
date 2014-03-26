"""
A viewer for astrodendro that looks like an integrated l,b and l,v thing.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

class IntegratedViewer(object):
    def __init__(self, dendrogram, hub):

        if dendrogram.data.ndim != 3:
            raise ValueError(
                "Only 3-dimensional arrays are supported")

        self.hub = hub # good
        self.hub.add_callback(self.update_selection) # good

        self.datacube = dendrogram.data
        self.dendrogram = dendrogram # good

        self.selected_contours = {} # selection_id -> (contour, contour) tuple

        self.fig = plt.figure(figsize=(10, 4.4))

        self.ax_lb = self.fig.add_axes([0.1, 0.05, 0.8, 0.4])
        self.ax_lv = self.fig.add_axes([0.1, 0.5, 0.8, 0.5])

        self.array_lb = np.nansum(self.datacube, axis=0)
        self.array_lv = np.nansum(self.datacube, axis=1)

        self._draw_plot()
        self.hub.add_callback(self.update_selection)

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

    def update_selection(self, selection_id):
        """Highlight seleted structures"""

        self.remove_contour(selection_id)
        
        if selection_id in self.hub.selections:

            struct = self.hub.selections[selection_id]

            if len(struct) != 1:
                raise NotImplemented(
                    "Multiple structures per selection not supported")

            struct = struct[0]

            if struct is not None:
                mask = struct.get_mask(subtree=True)
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