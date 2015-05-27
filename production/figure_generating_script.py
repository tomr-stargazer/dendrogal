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

from dendrogal.integrated_viewer import IntegratedViewer

from dendrogal.production.cloud_extractor_q1 import first_quad_cloud_catalog, compile_firstquad_catalog, export_firstquad_catalog
from dendrogal.production.cloud_extractor_q1 import d as quad1_d

from dendrogal.production.cloud_extractor_q2 import export_secondquad_catalog
from dendrogal.production.cloud_extractor_q2 import d as quad2_d

from dendrogal.production.cloud_extractor_q3 import export_thirdquad_catalog
from dendrogal.production.cloud_extractor_q3 import d as quad3_d

from dendrogal.production.cloud_extractor_q4 import fourth_quad_cloud_catalog, compile_fourthquad_catalog, export_fourthquad_catalog
from dendrogal.production.cloud_extractor_q4 import d as quad4_d

from dendrogal.production.cloud_extractor_carina import carina_cloud_catalog, compile_carina_catalog
from dendrogal.production.cloud_extractor_carina import d as carina_d

from dendrogal.comparison_to_other_catalogs import plot_dame_ellipses_on_integrated_viewer, plot_dame_ellipses_on_imf
from dendrogal.production.remove_degenerate_structures import reduce_catalog, selection_from_catalog
from dendrogal.production.map_figures import make_quadrant_lbv_map, make_quadrant_topdown_map, make_lv_map_new
from dendrogal.production.plot_catalog_measurements import plot_cmf, plot_size_linewidth_fit
from dendrogal.production.multipanel_catalog_measurements import multipanel_size_linewidth, multipanel_cmf

from dendrogal.production.map_dendrogram_thumbnail_figure import make_thumbnail_dendro_figure


output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def first_quadrant_figures(args=None, save=True):
    """ Creates the figures for the first quadrant and saves them. """

    cloud_catalog = compile_firstquad_catalog(first_quad_cloud_catalog())
    d = quad1_d

    imf_dame_a = make_lv_map_new(cloud_catalog, d)
    imf_dame_a.ax.set_xlim(135, 496)
    imf_dame_a.ax.set_ylim(138, 380)
    plot_dame_ellipses_on_imf(imf_dame_a)
    imf_dame_a.fig.savefig(output_path+"quad1_comparison_a.pdf", bbox_inches='tight')

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

    cmf_multipanel = multipanel_cmf()
    cmf_multipanel.savefig(output_path+"multipanel_cmf.pdf", bbox_inches='tight')

    size_linewidth_multipanel = multipanel_size_linewidth()
    size_linewidth_multipanel.savefig(output_path+'multipanel_size_linewidth.pdf', bbox_inches='tight')

    size_linewidth_allquad = plot_size_linewidth_fit(catalog)[0]
    plt.xlim(6, 300)
    plt.ylim(0.55, 13)
    size_linewidth_allquad.savefig(output_path+"allquads_size_linewidth.pdf", bbox_inches='tight')   
     

def make_thumbnail_figures():

    quad1_cloud_idx = 2772
    quad1_cat = export_firstquad_catalog()

    quad1_fig = make_thumbnail_dendro_figure(quad1_d, quad1_cat, quad1_cloud_idx)
    quad1_fig.savefig(output_path+"quad1_thumbnail.pdf", bbox_inches='tight')


    quad2_cloud_idx = 45
    quad2_cat = export_secondquad_catalog()

    quad2_fig = make_thumbnail_dendro_figure(quad2_d, quad2_cat, quad2_cloud_idx)
    quad2_fig.savefig(output_path+"quad2_thumbnail.pdf", bbox_inches='tight')

    # pseudo 3rd quadrant
    # quad2b_cloud_idx = 1506
    # quad3_fig = make_thumbnail_dendro_figure(quad2_d, quad2_cat, quad2b_cloud_idx)
    # quad3_fig.savefig(output_path+"quad3_thumbnail.pdf", bbox_inches='tight')

    quad3_cloud_idx = 726
    quad3_cat = export_thirdquad_catalog()

    quad3_fig = make_thumbnail_dendro_figure(quad3_d, quad3_cat, quad3_cloud_idx)
    # quad3_fig.ax_lb.set_ylim(0+2, 45-2)
    quad3_fig.savefig(output_path+"quad3_thumbnail.pdf", bbox_inches='tight')


    quad4_cloud_idx = 947
    quad4_cat = export_fourthquad_catalog()

    quad4_fig = make_thumbnail_dendro_figure(quad4_d, quad4_cat, quad4_cloud_idx)
    quad4_fig.savefig(output_path+"quad4_thumbnail.pdf", bbox_inches='tight')



