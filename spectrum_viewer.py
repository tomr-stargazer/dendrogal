"""
A spectrum viewer for astrodendro.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from wcsaxes import WCSAxes

class SpectrumViewer(object):
	def __init__(self, dendrogram, hub):

		self.hub = hub

        self.datacube = dendrogram.data
        self.dendrogram = dendrogram

        self.fig = plt.figure()

        ax = WCSAxes(self.fig, [0.1, 0.1, 0.8, 0.8] wcs=self.dendrogram.wcs, slices=(0,0,'x'))

        self.ax = self.fig.add_axes(ax)
