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

from astrodendro_analysis.integrated_viewer import IntegratedViewer

from astrodendro_analysis.production.new_cloud_extractor_q1 import first_quad_cloud_catalog, compile_firstquad_catalog
from astrodendro_analysis.production.new_cloud_extractor_q1 import d as quad1_d

from astrodendro_analysis.production.new_cloud_extractor_q2 import export_secondquad_catalog
from astrodendro_analysis.production.new_cloud_extractor_q2 import d as quad2_d

from astrodendro_analysis.production.new_cloud_extractor_q3 import export_thirdquad_catalog
from astrodendro_analysis.production.new_cloud_extractor_q3 import d as quad3_d

from astrodendro_analysis.production.new_cloud_extractor_q4 import fourth_quad_cloud_catalog, compile_fourthquad_catalog
from astrodendro_analysis.production.new_cloud_extractor_q4 import d as quad4_d

from astrodendro_analysis.production.new_cloud_extractor_carina import carina_cloud_catalog, compile_carina_catalog
from astrodendro_analysis.production.new_cloud_extractor_carina import d as carina_d

from astrodendro_analysis.comparison_to_other_catalogs import plot_dame_ellipses_on_integrated_viewer, plot_dame_ellipses_on_imf
from astrodendro_analysis.production.remove_degenerate_structures import reduce_catalog, selection_from_catalog
from astrodendro_analysis.production.map_figures import make_quadrant_lbv_map, make_quadrant_topdown_map, make_lv_map_new
from astrodendro_analysis.production.plot_catalog_measurements import plot_cmf, plot_size_linewidth_fit


output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def first_quadrant_figures(args=None, save=True):
    """ Creates the figures for the first quadrant and saves them. """

    cloud_catalog = compile_firstquad_catalog(first_quad_cloud_catalog())
    d = quad1_d

    # cloud_selection = selection_from_catalog(d, cloud_catalog)

    imf_dame_a = make_lv_map_new(cloud_catalog, d)
    imf_dame_a.ax.set_xlim(135, 496)
    imf_dame_a.ax.set_ylim(138, 380)
    plot_dame_ellipses_on_imf(imf_dame_a)
    imf_dame_a.fig.savefig(output_path+"quad1_comparison_a.pdf", bbox_inches='tight')
    # # make the first one
    # iv_dame_a = make_quadrant_lbv_map(cloud_catalog, d)
    # iv_dame_a.ax_lv.set_xlim(115, 496.5)
    # iv_dame_a.ax_lv.set_ylim(120, 413)
    # iv_dame_a.fig.savefig(output_path+"quad1_comparison_a.pdf", bbox_inches='tight')

    # iv_dame_b = make_quadrant_lbv_map(cloud_catalog, d, contour_select=False)
    # plot_dame_ellipses_on_integrated_viewer(iv_dame_b)
    # iv_dame_b.ax_lv.set_xlim(115, 496.5)
    # iv_dame_b.ax_lv.set_ylim(120, 413)

    # iv_dame_b.fig.savefig(output_path+"quad1_comparison_b.pdf", bbox_inches='tight')

    imf_firstquad = make_lv_map_new(cloud_catalog, d)
    # imf_firstquad.fig.set_figheight(5.7)
    # imf_firstquad.fig.set_figwidth(10)
    imf_firstquad.fig.canvas.draw()
    imf_firstquad.fig.savefig(output_path+'quad1_map.pdf', bbox_inches='tight')

    topdown_1q = make_quadrant_topdown_map(cloud_catalog, figsize=(6,8))
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

    cloud_catalog = export_secondquad_catalog()
    d = quad2_d

    # cloud_selection = selection_from_catalog(d, cloud_catalog)

    # iv_secondquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal', aspect=2, vscale=1.3)
    # iv_secondquad.fig.set_figheight(5.7)
    # iv_secondquad.fig.set_figwidth(10)
    # iv_secondquad.fig.canvas.draw()
    # iv_secondquad.fig.savefig(output_path+'quad2_map.pdf', bbox_inches='tight')

    imf_secondquad = make_lv_map_new(cloud_catalog, d, integration_limits=(-2,2))
    imf_secondquad.ax.coords['glon'].set_ticks(spacing=10*u.deg, color='white', exclude_overlapping=True)    
    imf_secondquad.ax.set_xlim(30,1100)
    imf_secondquad.ax.set_ylim(30,145)
    # imf_secondquad.fig.set_figheight(5.7)
    # imf_secondquad.fig.set_figwidth(10)
    imf_secondquad.fig.canvas.draw()
    imf_secondquad.fig.savefig(output_path+'quad2_map.pdf', bbox_inches='tight')

    topdown_2q = make_quadrant_topdown_map(cloud_catalog, loc='lower right', figsize=(6,8))
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

    cloud_catalog = export_thirdquad_catalog()
    d = quad3_d

    # iv_thirdquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal', aspect=2, vscale=1.3)
    # iv_thirdquad.fig.set_figheight(5.7)
    # iv_thirdquad.fig.set_figwidth(10)
    # iv_thirdquad.fig.canvas.draw()
    # iv_thirdquad.fig.savefig(output_path+'quad3_map.pdf', bbox_inches='tight')

    imf_thirdquad = make_lv_map_new(cloud_catalog, d, integration_limits=(-1.9,3))
    imf_thirdquad.ax.coords['glon'].set_ticks(spacing=5*u.deg, color='white', exclude_overlapping=True)
    imf_thirdquad.ax.set_xlim(135, 744)
    imf_thirdquad.fig.canvas.draw()
    imf_thirdquad.fig.savefig(output_path+'quad3_map.pdf', bbox_inches='tight')

    topdown_3q = make_quadrant_topdown_map(cloud_catalog, loc='lower right', figsize=(6,8))
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


    cloud_catalog = compile_fourthquad_catalog(fourth_quad_cloud_catalog())
    d = quad4_d


    # # if args is None, the dendrogram will be computed fresh
    # if args is None:
    #     d, catalog, header, metadata = fourth_quad_dendrogram()
    #     args = (d, catalog, header, metadata)
    # else:
    #     d, catalog, header, metadata = args

    # cloud_catalog = export_fourthquad_catalog(args=args)

    # cloud_selection = selection_from_catalog(d, cloud_catalog)

    # iv_fourthquad = make_quadrant_lbv_map(cloud_catalog, d, alignment='horizontal', aspect=1/2, vscale=1.3)
    # iv_fourthquad.fig.set_figheight(8)
    # iv_fourthquad.fig.set_figwidth(8)
    # iv_fourthquad.fig.canvas.draw()
    # iv_fourthquad.fig.savefig(output_path+'quad4_map.pdf', bbox_inches='tight')

    imf_fourthquad = make_lv_map_new(cloud_catalog, d)
    imf_fourthquad.ax.set_ylim(20,245) # found via trial & error
    imf_fourthquad.fig.canvas.draw()
    imf_fourthquad.fig.savefig(output_path+'quad4_map.pdf', bbox_inches='tight')

    topdown_4q = make_quadrant_topdown_map(cloud_catalog, loc='lower right', figsize=(6,8))
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


def carina_figures(args=None, save=True):
    """ Creates the figures for the fourth quadrant and saves them. """


    cloud_catalog = compile_carina_catalog(carina_cloud_catalog())
    d = carina_d

    imf_carina = make_lv_map_new(cloud_catalog, d)
    imf_carina.ax.set_ylim(70,200)
    imf_carina.fig.canvas.draw()
    imf_carina.fig.savefig(output_path+'carina_map.pdf', bbox_inches='tight')

    topdown_carina = make_quadrant_topdown_map(cloud_catalog, loc='lower right', figsize=(6,8))
    topdown_carina.axes[0].set_xlim(-13.1, 5.1)
    topdown_carina.axes[0].set_ylim(17.6, -0.05)
    topdown_carina.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_carina.axes[0].set_ylabel("Solar-centric $x$ (kpc)")

    topdown_carina.savefig(output_path+'carina_topdown.pdf', bbox_inches='tight')

    cmf_carina = plot_cmf(cloud_catalog)[0]
    cmf_carina.savefig(output_path+"carina_cmf.pdf", bbox_inches='tight')

    size_linewidth_carina = plot_size_linewidth_fit(cloud_catalog)[0]
    plt.xlim(3, 115)
    plt.ylim(0.55, 29)
    size_linewidth_carina.savefig(output_path+"carina_size_linewidth.pdf", bbox_inches='tight')    

def combined_galaxy_figures(catalog):

    # catalog = extract_and_combine_catalogs(args)

    topdown_all = make_quadrant_topdown_map(catalog, loc='lower center', figsize=(9,9))
    topdown_all.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_all.axes[0].set_ylabel("Solar-centric $x$ (kpc)")
    topdown_all.axes[0].set_xlim(-13, 13)
    topdown_all.axes[0].set_ylim(21.9, -6.5)

    topdown_all.savefig(output_path+"allquads_topdown.pdf", bbox_inches='tight')

    cmf_allquad = plot_cmf(catalog)[0]
    cmf_allquad.savefig(output_path+"allquads_cmf.pdf", bbox_inches='tight')

    size_linewidth_allquad = plot_size_linewidth_fit(catalog)[0]
    plt.xlim(6, 300)
    plt.ylim(0.55, 13)
    size_linewidth_allquad.savefig(output_path+"allquads_size_linewidth.pdf", bbox_inches='tight')   
     





