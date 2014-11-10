"""
A spectrum viewer for astrodendro.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from wcsaxes import WCSAxes

class SpectrumViewer(object):
    def __init__(self, dendrogram, hub, l, b):

        self.hub = hub

        # coordinates to plot
        self.l = l
        self.b = b 

        self.datacube = dendrogram.data
        self.dendrogram = dendrogram

        self.fig = plt.figure(figsize=(11.5,3.7))
        # ax = WCSAxes(self.fig, [0.1, 0.1, 0.8, 0.8], wcs=self.dendrogram.wcs, slices=(0,0,'x'))
        # ax = plt.subplot(111)
        self.ax = self.fig.add_axes([0.1, 0.15, 0.8, 0.75])
        # self.ax = self.fig.add_axes(ax)

        self._draw_plot()

    def _draw_plot(self):

        spectrum = self.datacube[:, self.b, self.l]

        pixels = np.vstack([np.ones(len(spectrum))*self.l, np.ones(len(spectrum))*self.b, np.arange(len(spectrum))]).T

        coordinates = self.dendrogram.wcs.all_pix2world(pixels, 0)
        velocities = coordinates[:,2]

        self.ax.plot(velocities, spectrum, 'k')
        self.ax.set_xlim(min(velocities), max(velocities))
        self.ax.set_xlabel(r"v$_{LSR}$ (km s$^{-1}$)")
        self.ax.set_ylabel(r"T$_{MB}$ (K)")

        self.fig.canvas.draw()