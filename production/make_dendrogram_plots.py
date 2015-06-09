"""
This makes dendrogram plots which could be applied to the quadrant catalogs.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from dendrogal.integrated_viewer import IntegratedViewer
from dendrogal.production.remove_degenerate_structures import reduce_catalog, selection_from_catalog
from dendrogal.hurt_image_overlay_demo import underlay_hurt_galaxy
from dendrogal.production.integrated_map_figure import IntegratedMapFigure

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'

def quadrant_dendrogram_with_highlights(d, catalog):
    """ Plots a quadrant's dendrogram, with the catalog highlighted. """

    fig = plt.figure()

    ax = fig.add_subplot(111)

    p = d.plotter()

    # get the whole tree up on there
    p.plot_tree(ax, color='black')

    # now do each thing individually
    for item in catalog:
        idx = item['_idx']

        p.plot_tree(ax, structure=[d[idx]], color=colorbrewer_red, lw=2)

    fig.canvas.draw()

    return fig
