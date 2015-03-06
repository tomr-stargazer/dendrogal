"""
This script generates the figures that end up in the PDF paper document.

It does this slowly.

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

from ..integrated_viewer import IntegratedViewer

# first quadrant figures
from .first_quadrant_cloud_extraction import first_quad_dendrogram, export_firstquad_catalog
from ..comparison_to_other_catalogs import plot_dame_ellipses_on_integrated_viewer
from .remove_degenerate_structures import reduce_catalog, selection_from_catalog
from .map_figures import make_quadrant_lbv_map, make_quadrant_topdown_map
from .plot_catalog_measurements import plot_cmf, plot_size_linewidth


output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def first_quadrant_figures(args=None, save=True):
    """ Creates the figures for the first quadrant and saves them. """

    # if args is None, the dendrogram will be computed fresh
    if args is None:
        d, catalog, header, metadata = first_quad_dendrogram()
        args = (d, catalog, header, metadata)
    else:
        d, catalog, header, metadata = args

    cloud_catalog = export_firstquad_catalog(args=args)

    cloud_selection = selection_from_catalog(d, cloud_catalog)

    # make the first one
    iv_dame_a = make_quadrant_lbv_map(cloud_catalog, d)
    iv_dame_a.ax_lv.set_xlim(115, 496.5)
    iv_dame_a.ax_lv.set_ylim(120, 413)
    iv_dame_a.fig.savefig(output_path+"quad1_comparison_a.pdf", bbox_inches='tight')

    iv_dame_b = make_quadrant_lbv_map(cloud_catalog, d, contour_select=False)
    plot_dame_ellipses_on_integrated_viewer(iv_dame_b)
    iv_dame_b.ax_lv.set_xlim(115, 496.5)
    iv_dame_b.ax_lv.set_ylim(120, 413)

    iv_dame_b.fig.savefig(output_path+"quad1_comparison_b.pdf", bbox_inches='tight')

    iv_firstquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal')
    iv_firstquad.fig.set_figheight(5.7)
    iv_firstquad.fig.set_figwidth(10)
    iv_firstquad.fig.canvas.draw()
    iv_firstquad.fig.savefig(output_path+'quad1_map.pdf')

    topdown_1q = make_quadrant_topdown_map(cloud_catalog)
    topdown_1q.set_figwidth(6)
    topdown_1q.set_figheight(8)
    topdown_1q.axes[0].set_xlim(0, 15)
    topdown_1q.axes[0].set_ylim(17, -4.5)
    topdown_1q.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_1q.axes[0].set_ylabel("Solar-centric $x$ (kpc)")

    topdown_1q.savefig(output_path+'quad1_topdown.pdf')

    cmf_1q = plot_cmf(cloud_catalog)[0]
    cmf_1q.savefig(output_path+"quad1_cmf.pdf")

    size_linewidth_1q = plot_size_linewidth(cloud_catalog)
    plt.xlim(7, 115)
    plt.ylim(0.55, 29)
    size_linewidth_1q.savefig(output_path+"quad1_size_linewith.pdf")

def second_quadrant_figures(args=None, save=True):
    """ Creates the figures for the second quadrant and saves them. """

    # if args is None, the dendrogram will be computed fresh
    if args is None:
        d, catalog, header, metadata = second_quad_dendrogram()
        args = (d, catalog, header, metadata)
    else:
        d, catalog, header, metadata = args

    cloud_catalog = export_secondquad_catalog(args=args)

    cloud_selection = selection_from_catalog(d, cloud_catalog)

    iv_firstquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal')
    iv_firstquad.fig.set_figheight(5.7)
    iv_firstquad.fig.set_figwidth(10)
    iv_firstquad.fig.canvas.draw()
    iv_firstquad.fig.savefig(output_path+'quad2_map.pdf')

    topdown_1q = make_quadrant_topdown_map(cloud_catalog)
    topdown_1q.set_figwidth(6)
    topdown_1q.set_figheight(8)
    topdown_1q.axes[0].set_xlim(0, 15)
    topdown_1q.axes[0].set_ylim(17, -4.5)
    topdown_1q.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_1q.axes[0].set_ylabel("Solar-centric $x$ (kpc)")

    topdown_1q.savefig(output_path+'quad2_topdown.pdf')

    cmf_1q = plot_cmf(cloud_catalog)[0]
    cmf_1q.savefig(output_path+"quad2_cmf.pdf")

    size_linewidth_1q = plot_size_linewidth(cloud_catalog)
    plt.xlim(7, 115)
    plt.ylim(0.55, 29)
    size_linewidth_1q.savefig(output_path+"quad2_size_linewith.pdf")

