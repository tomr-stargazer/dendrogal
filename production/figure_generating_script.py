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

from .first_quadrant_cloud_extraction import first_quad_dendrogram, export_firstquad_catalog
from .second_quadrant_cloud_extraction import second_quad_dendrogram, export_secondquad_catalog
from .third_quadrant_cloud_extraction import third_quad_dendrogram, export_thirdquad_catalog
from .fourth_quadrant_cloud_extraction import fourth_quad_dendrogram, export_fourthquad_catalog
from ..comparison_to_other_catalogs import plot_dame_ellipses_on_integrated_viewer
from .remove_degenerate_structures import reduce_catalog, selection_from_catalog
from .map_figures import make_quadrant_lbv_map, make_quadrant_topdown_map, make_lv_map_new
from .plot_catalog_measurements import plot_cmf, plot_size_linewidth_fit


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

    # iv_firstquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal')
    # iv_firstquad.fig.set_figheight(5.7)
    # iv_firstquad.fig.set_figwidth(10)
    # iv_firstquad.fig.canvas.draw()
    # iv_firstquad.fig.savefig(output_path+'quad1_map.pdf', bbox_inches='tight')

    imf_firstquad = make_lv_map_new(cloud_catalog, d)
    # imf_firstquad.fig.set_figheight(5.7)
    # imf_firstquad.fig.set_figwidth(10)
    imf_firstquad.fig.canvas.draw()
    imf_firstquad.fig.savefig(output_path+'quad1_map.pdf', bbox_inches='tight')

    topdown_1q = make_quadrant_topdown_map(cloud_catalog)
    topdown_1q.set_figwidth(6)
    topdown_1q.set_figheight(8)
    topdown_1q.axes[0].set_xlim(0, 15)
    topdown_1q.axes[0].set_ylim(17, -4.5)
    topdown_1q.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_1q.axes[0].set_ylabel("Solar-centric $x$ (kpc)")

    topdown_1q.savefig(output_path+'quad1_topdown.pdf', bbox_inches='tight')

    cmf_1q = plot_cmf(cloud_catalog)[0]
    cmf_1q.savefig(output_path+"quad1_cmf.pdf", bbox_inches='tight')

    size_linewidth_1q = plot_size_linewidth_fit(cloud_catalog)[0]
    plt.xlim(3, 115)
    plt.ylim(0.55, 29)
    size_linewidth_1q.savefig(output_path+"quad1_size_linewidth.pdf", bbox_inches='tight')

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

    # iv_secondquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal', aspect=2, vscale=1.3)
    # iv_secondquad.fig.set_figheight(5.7)
    # iv_secondquad.fig.set_figwidth(10)
    # iv_secondquad.fig.canvas.draw()
    # iv_secondquad.fig.savefig(output_path+'quad2_map.pdf', bbox_inches='tight')

    imf_secondquad = make_lv_map_new(cloud_catalog, d, integration_range=(-2,2))
    # imf_secondquad.fig.set_figheight(5.7)
    # imf_secondquad.fig.set_figwidth(10)
    imf_secondquad.fig.canvas.draw()
    imf_secondquad.fig.savefig(output_path+'quad2_map.pdf', bbox_inches='tight')

    topdown_2q = make_quadrant_topdown_map(cloud_catalog, loc='lower right')
    topdown_2q.set_figwidth(6)
    topdown_2q.set_figheight(8)
    topdown_2q.axes[0].set_xlim(-3, 12)
    topdown_2q.axes[0].set_ylim(8.3, -7.6)
    topdown_2q.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_2q.axes[0].set_ylabel("Solar-centric $x$ (kpc)")

    topdown_2q.savefig(output_path+'quad2_topdown.pdf', bbox_inches='tight')

    cmf_2q = plot_cmf(cloud_catalog, min_mass=1e4)[0]
    cmf_2q.savefig(output_path+"quad2_cmf.pdf", bbox_inches='tight')

    size_linewidth_2q = plot_size_linewidth_fit(cloud_catalog)[0]
    plt.xlim(3, 115)
    plt.ylim(0.55, 29)
    size_linewidth_2q.savefig(output_path+"quad2_size_linewidth.pdf", bbox_inches='tight')

def third_quadrant_figures(args=None, save=True):
    """ Creates the figures for the third quadrant and saves them. """

    # if args is None, the dendrogram will be computed fresh
    if args is None:
        d, catalog, header, metadata = third_quad_dendrogram()
        args = (d, catalog, header, metadata)
    else:
        d, catalog, header, metadata = args

    cloud_catalog = export_thirdquad_catalog(args=args)

    cloud_selection = selection_from_catalog(d, cloud_catalog)

    # iv_thirdquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal', aspect=2, vscale=1.3)
    # iv_thirdquad.fig.set_figheight(5.7)
    # iv_thirdquad.fig.set_figwidth(10)
    # iv_thirdquad.fig.canvas.draw()
    # iv_thirdquad.fig.savefig(output_path+'quad3_map.pdf', bbox_inches='tight')

    imf_thirdquad = make_lv_map_new(cloud_catalog, d, integration_range=(-2,2))
    imf_thirdquad.fig.canvas.draw()
    imf_thirdquad.fig.savefig(output_path+'quad3_map.pdf', bbox_inches='tight')

    topdown_3q = make_quadrant_topdown_map(cloud_catalog, loc='lower right')
    topdown_3q.set_figwidth(6)
    topdown_3q.set_figheight(8)
    topdown_3q.axes[0].set_xlim(-11.4, 1.6)
    topdown_3q.axes[0].set_ylim(5.3, -6.7)
    topdown_3q.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_3q.axes[0].set_ylabel("Solar-centric $x$ (kpc)")

    topdown_3q.savefig(output_path+'quad3_topdown.pdf', bbox_inches='tight')

    cmf_3q = plot_cmf(cloud_catalog, min_mass=1e4)[0]
    cmf_3q.savefig(output_path+"quad3_cmf.pdf", bbox_inches='tight')

    size_linewidth_3q = plot_size_linewidth_fit(cloud_catalog)[0]
    plt.xlim(3, 115)
    plt.ylim(0.55, 29)
    size_linewidth_3q.savefig(output_path+"quad3_size_linewidth.pdf", bbox_inches='tight')    

def fourth_quadrant_figures(args=None, save=True):
    """ Creates the figures for the fourth quadrant and saves them. """

    # if args is None, the dendrogram will be computed fresh
    if args is None:
        d, catalog, header, metadata = fourth_quad_dendrogram()
        args = (d, catalog, header, metadata)
    else:
        d, catalog, header, metadata = args

    cloud_catalog = export_fourthquad_catalog(args=args)

    cloud_selection = selection_from_catalog(d, cloud_catalog)

    # iv_fourthquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal', aspect=1/2, vscale=1.3)
    # iv_fourthquad.fig.set_figheight(8)
    # iv_fourthquad.fig.set_figwidth(8)
    # iv_fourthquad.fig.canvas.draw()
    # iv_fourthquad.fig.savefig(output_path+'quad4_map.pdf', bbox_inches='tight')

    imf_fourthquad = make_lv_map_new(cloud_catalog, d)
    imf_fourthquad.fig.canvas.draw()
    imf_fourthquad.fig.savefig(output_path+'quad4_map.pdf', bbox_inches='tight')

    topdown_4q = make_quadrant_topdown_map(cloud_catalog, loc='lower right')
    topdown_4q.set_figwidth(6)
    topdown_4q.set_figheight(8)
    topdown_4q.axes[0].set_xlim(-13.1, 5.1)
    topdown_4q.axes[0].set_ylim(17.6, -0.05)
    topdown_4q.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_4q.axes[0].set_ylabel("Solar-centric $x$ (kpc)")

    topdown_4q.savefig(output_path+'quad4_topdown.pdf', bbox_inches='tight')

    cmf_4q = plot_cmf(cloud_catalog)[0]
    cmf_4q.savefig(output_path+"quad4_cmf.pdf", bbox_inches='tight')

    size_linewidth_4q = plot_size_linewidth_fit(cloud_catalog)[0]
    plt.xlim(3, 115)
    plt.ylim(0.55, 29)
    size_linewidth_4q.savefig(output_path+"quad4_size_linewidth.pdf", bbox_inches='tight')    

def combined_galaxy_figures(catalog):

    # catalog = extract_and_combine_catalogs(args)

    topdown_all = make_quadrant_topdown_map(catalog, loc='lower center')
    topdown_all.set_figwidth(9)
    topdown_all.set_figheight(9)
    topdown_all.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_all.axes[0].set_ylabel("Solar-centric $x$ (kpc)")
    topdown_all.axes[0].set_xlim(-13, 13)
    topdown_all.axes[0].set_ylim(21.9, -6.5)

    topdown_all.savefig(output_path+"allquads_topdown.pdf", bbox_inches='tight')

    cmf_allquad = plot_cmf(catalog)[0]
    cmf_allquad.savefig(output_path+"allquads_cmf.pdf", bbox_inches='tight')

    size_linewidth_allquad = plot_size_linewidth_fit(catalog)[0]
    plt.xlim(3, 115)
    plt.ylim(0.55, 29)
    size_linewidth_allquad.savefig(output_path+"allquads_size_linewidth.pdf", bbox_inches='tight')   
     





