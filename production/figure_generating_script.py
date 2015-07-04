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
from dendrogal.production.map_figures import make_quadrant_lbv_map, make_quadrant_topdown_map, make_lv_map_new, draw_solar_and_tangent_circles
from dendrogal.production.plot_catalog_measurements import plot_cmf, plot_size_linewidth_fit
from dendrogal.production.multipanel_catalog_measurements import multipanel_size_linewidth, multipanel_cmf

from dendrogal.production.map_dendrogram_thumbnail_figure import make_thumbnail_dendro_figure


output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def first_quadrant_figures(args=None, save=True):
    """ Creates the figures for the first quadrant and saves them. """

    cloud_catalog = compile_firstquad_catalog(first_quad_cloud_catalog())
    d = quad1_d

    imf_firstquad = make_lv_map_new(cloud_catalog, d, integration_limits=(-2, 2), ellipse_thickness_inner=1.2, ellipse_thickness_outer=2.4)
    imf_firstquad.ax.set_ylim(15, 392)    
    plot_dame_ellipses_on_imf(imf_firstquad)    
    imf_firstquad.fig.canvas.draw()
    imf_firstquad.fig.savefig(output_path+'quad1_map.pdf', bbox_inches='tight')


def second_quadrant_figures(args=None, save=True):
    """ Creates the figures for the second quadrant and saves them. """

    cloud_catalog = export_secondquad_catalog()
    d = quad2_d

    imf_secondquad = make_lv_map_new(cloud_catalog, d, integration_limits=(-2,2), ellipse_thickness_inner=0.6, ellipse_thickness_outer=1.4)
    imf_secondquad.ax.coords['glon'].set_ticks(spacing=10*u.deg, color='white', exclude_overlapping=True)    
    imf_secondquad.ax.set_xlim(30,1100)
    imf_secondquad.ax.set_ylim(30,145)
    imf_secondquad.fig.canvas.draw()
    imf_secondquad.fig.savefig(output_path+'quad2_map.pdf', bbox_inches='tight')


def third_quadrant_figures(args=None, save=True):
    """ Creates the figures for the third quadrant and saves them. """

    cloud_catalog = export_thirdquad_catalog()
    d = quad3_d

    # these bizarre "integration limits" have to do with a bug in the WCS for the third quadrant.
    imf_thirdquad = make_lv_map_new(cloud_catalog, d, integration_limits=(-1.9,3), ellipse_thickness_inner=0.6, ellipse_thickness_outer=1.4)
    imf_thirdquad.ax.coords['glon'].set_ticks(spacing=5*u.deg, color='white', exclude_overlapping=True)
    imf_thirdquad.ax.set_xlim(135, 744)
    imf_thirdquad.fig.canvas.draw()
    imf_thirdquad.fig.savefig(output_path+'quad3_map.pdf', bbox_inches='tight')


def fourth_quadrant_figures(args=None, save=True):
    """ Creates the figures for the fourth quadrant and saves them. """


    cloud_catalog = compile_fourthquad_catalog(fourth_quad_cloud_catalog())
    d = quad4_d

    imf_fourthquad = make_lv_map_new(cloud_catalog, d, integration_limits=(-2, 2))
    imf_fourthquad.ax.set_ylim(45,216) # found via trial & error
    imf_fourthquad.fig.canvas.draw()
    imf_fourthquad.fig.savefig(output_path+'quad4_map.pdf', bbox_inches='tight')


def carina_figures(args=None, save=True):
    """ Creates the figures for the fourth quadrant and saves them. """


    cloud_catalog = compile_carina_catalog(carina_cloud_catalog())
    d = carina_d

    imf_carina = make_lv_map_new(cloud_catalog, d, integration_limits=(-2, 2), ellipse_thickness_inner=1.5, ellipse_thickness_outer=2.5)
    imf_carina.ax.set_ylim(88,170)
    imf_carina.ax.set_xlim(-0.5,160.1)
    imf_carina.fig.canvas.draw()
    imf_carina.fig.savefig(output_path+'carina_map.pdf', bbox_inches='tight')


def combined_galaxy_figures(catalog):

    # catalog = extract_and_combine_catalogs(args)

    topdown_all = make_quadrant_topdown_map(catalog, loc='lower center', figsize=(9,9))
    draw_solar_and_tangent_circles(topdown_all)    
    plt.text(9, 17, 'IQ', fontsize=20, family='serif', color='w')
    plt.text(9, -4, 'IIQ', fontsize=20, family='serif', color='w')
    plt.text(-11.5, -4, 'IIIQ', fontsize=20, family='serif', color='w')
    plt.text(-11.5, 17, 'IVQ', fontsize=20, family='serif', color='w')

    topdown_all.axes[0].set_xlabel("Solar-centric $y$ (kpc)")
    topdown_all.axes[0].set_ylabel("Solar-centric $x$ (kpc)")
    topdown_all.axes[0].set_xlim(-13, 13)
    topdown_all.axes[0].set_ylim(21.9, -6.5)

    topdown_all.savefig(output_path+"allquads_topdown.pdf", bbox_inches='tight')

    cmf_multipanel = multipanel_cmf()
    cmf_multipanel.savefig(output_path+"multipanel_cmf.pdf", bbox_inches='tight')

    size_linewidth_multipanel = multipanel_size_linewidth()
    size_linewidth_multipanel.savefig(output_path+'multipanel_size_linewidth.pdf', bbox_inches='tight')
     

def make_thumbnail_figures():

    quad1_cloud_idx = 2868
    quad1_cat = export_firstquad_catalog()

    quad1_fig = make_thumbnail_dendro_figure(quad1_d, quad1_cat, quad1_cloud_idx)
    quad1_fig.savefig(output_path+"quad1_thumbnail.pdf", bbox_inches='tight')

    print quad1_cat['mass'][quad1_cat['_idx']==quad1_cloud_idx]
    print quad1_cat['error_mass_plus'][quad1_cat['_idx']==quad1_cloud_idx]
    print quad1_cat['error_mass_minus'][quad1_cat['_idx']==quad1_cloud_idx]

    print quad1_cat['distance'][quad1_cat['_idx']==quad1_cloud_idx]
    print quad1_cat['error_distance_plus'][quad1_cat['_idx']==quad1_cloud_idx]
    print quad1_cat['error_distance_minus'][quad1_cat['_idx']==quad1_cloud_idx]

    print ""

    quad2_cloud_idx = 66
    quad2_cat = export_secondquad_catalog()

    quad2_fig = make_thumbnail_dendro_figure(quad2_d, quad2_cat, quad2_cloud_idx)
    quad2_fig.savefig(output_path+"quad2_thumbnail.pdf", bbox_inches='tight')

    print quad2_cat['mass'][quad2_cat['_idx']==quad2_cloud_idx]
    print quad2_cat['error_mass_plus'][quad2_cat['_idx']==quad2_cloud_idx]
    print quad2_cat['error_mass_minus'][quad2_cat['_idx']==quad2_cloud_idx]

    print quad2_cat['distance'][quad2_cat['_idx']==quad2_cloud_idx]
    print quad2_cat['error_distance_plus'][quad2_cat['_idx']==quad2_cloud_idx]
    print quad2_cat['error_distance_minus'][quad2_cat['_idx']==quad2_cloud_idx]

    print ""
    
    quad3_cloud_idx = 697
    quad3_cat = export_thirdquad_catalog()

    quad3_fig = make_thumbnail_dendro_figure(quad3_d, quad3_cat, quad3_cloud_idx, panel_width=3.5*u.deg, latitude_px_override=[22, 32])
    quad3_fig.savefig(output_path+"quad3_thumbnail.pdf", bbox_inches='tight')

    print quad3_cat['mass'][quad3_cat['_idx']==quad3_cloud_idx]
    print quad3_cat['error_mass_plus'][quad3_cat['_idx']==quad3_cloud_idx]
    print quad3_cat['error_mass_minus'][quad3_cat['_idx']==quad3_cloud_idx]

    print quad3_cat['distance'][quad3_cat['_idx']==quad3_cloud_idx]
    print quad3_cat['error_distance_plus'][quad3_cat['_idx']==quad3_cloud_idx]
    print quad3_cat['error_distance_minus'][quad3_cat['_idx']==quad3_cloud_idx]

    print ""


    quad3b_cloud_idx = 710
    quad3b_fig = make_thumbnail_dendro_figure(quad3_d, quad3_cat, quad3b_cloud_idx, latitude_px_override=[17, 34])
    quad3b_fig.savefig(output_path+"quad3b_thumbnail.pdf", bbox_inches='tight')

    quad4_cloud_idx = 484
    quad4_cat = export_fourthquad_catalog()

    quad4_fig = make_thumbnail_dendro_figure(quad4_d, quad4_cat, quad4_cloud_idx)
    quad4_fig.savefig(output_path+"quad4_thumbnail.pdf", bbox_inches='tight')

    print quad4_cat['mass'][quad4_cat['_idx']==quad4_cloud_idx]
    print quad4_cat['error_mass_plus'][quad4_cat['_idx']==quad4_cloud_idx]
    print quad4_cat['error_mass_minus'][quad4_cat['_idx']==quad4_cloud_idx]

    print quad4_cat['distance'][quad4_cat['_idx']==quad4_cloud_idx]
    print quad4_cat['error_distance_plus'][quad4_cat['_idx']==quad4_cloud_idx]
    print quad4_cat['error_distance_minus'][quad4_cat['_idx']==quad4_cloud_idx]


