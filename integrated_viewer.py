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

        self.selected_contour = {} # selection_id -> Contour

        self.fig = plt.figure(figsize=(14, 8))

        self.ax1 = self.fig.add_axes([0.1, 0.05, 0.8, 0.4]) # l, b
        self.ax2 = self.fig.add_axes([0.1, 0.5, 0.8, 0.5]) # l, v

        self.array_lb = np.nansum(self.datacube, axis=0)
        self.array_lv = np.nansum(self.datacube, axis=1)

        self._draw_plot()
        self.hub.add_callback(self.update_selection)

    def _draw_plot(self):
    	""" Create an image to plot things onto. """

        len_v, len_l, len_b = self.datacube.shape
        crop_factor = np.round(len_v/5)

        # Crop l,v plot because it normally goes awkwardly too high and low in v space
        self._clim_lb = (np.min(self.array_lb[~np.isnan(self.array_lb) & ~np.isinf(self.array_lb)]),
        	             np.max(self.array_lb[~np.isnan(self.array_lb) & ~np.isinf(self.array_lb)]))

        self._clim_lv = (np.min(self.array_lv[~np.isnan(self.array_lv) & ~np.isinf(self.array_lv)]),
        	             np.max(self.array_lv[~np.isnan(self.array_lv) & ~np.isinf(self.array_lv)]))

        self.image_lb = self.ax1.imshow(self.array_lb, origin='lower', 
            interpolation='nearest', vmin=self._clim_lb[0], 
            vmax=0.025*self._clim_lb[1], cmap=plt.cm.gray)

        self.image_lv = self.ax2.imshow(self.array_lv, origin='lower', 
            interpolation='nearest', vmin=self._clim_lv[0], 
            vmax=0.05*self._clim_lv[1], cmap=plt.cm.gray, aspect=2.5)
        self.ax2.set_ylim(crop_factor, len_v - crop_factor)

        self.fig.canvas.draw()

    def update_selection(self, selection_id):
        """Highlight seleted structures"""
        
        if selection_id in self.lines2d:
            if self.lines2d[selection_id] is not None:
                self.lines2d[selection_id].remove()
                del self.lines2d[selection_id]

        struct = self.hub.selections[selection_id][0]
        if struct is None:
            self.fig.canvas.draw()
            return
        selected_indices = [leaf.idx for leaf in struct.descendants + [struct]]

        self.lines2d[selection_id] = self.axes.plot(
            self.xdata[selected_indices], 
            self.ydata[selected_indices], 
            'o', color=self.hub.colors[selection_id], zorder=struct.height)[0]

        self.fig.canvas.draw()

    def select(self, structure, index):
        "Select a given structure and assign it to index-th selection"
        raise NotImplementedError()