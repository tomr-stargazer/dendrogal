"""
Uses a grid of cloud extraction parameters to optimize the first quadrant cloud extraction parameters.

"""

from __future__ import division

import numpy as np

from astrodendro_analysis import first_quadrant_cloud_extraction, map_figures
from astrodendro_analysis.comparison_to_other_catalogs import plot_dame_ellipses_on_imf

colorbrewer_red = '#e41a1c' 
colorbrewer_blue = '#377eb8'
colorbrewer_green = '#4daf4a'

min_npix_list = [20, 30, 40]
max_tree_list = [10, 30]
min_val_list = [0.18]
min_delta_list = [0.18]
min_tree_list = [2, 0]

print "Total upcoming iterations: {0}".format(len(max_tree_list)*len(min_val_list)*len(min_delta_list)*len(min_tree_list)*len(min_npix_list))

directory = 'grid_of_cloud_parameters/'

def run_metrics_on(cloud_catalog, dendrogram, 
                   min_val=None, min_delta=None, 
                   max_tree=None, min_tree=None, min_npix=None):

    # make and save some plots

    imf = map_figures.make_lv_map_new(cloud_catalog, dendrogram, ellipse_color=colorbrewer_blue)
    min_value_str = '{:.3f}'.format(min_val)
    min_delta_str = '{:.3f}'.format(min_delta)
    max_tree_str = str(max_tree)
    min_tree_str = str(min_tree)
    min_npix_str = str(min_npix)
    
    filename = "np{4}_mv{0}_md{1}_xt{2}_nt{3}.png".format(min_value_str, min_delta_str, max_tree_str, min_tree_str, min_npix_str)

    imf.fig.savefig(directory+filename, bbox_inches='tight')

    plot_dame_ellipses_on_imf(imf)

    imf.fig.savefig(directory+'dame_'+filename, bbox_inches='tight')



for min_val in min_val_list:

    for min_delta in min_delta_list:

        for min_npix in min_npix_list:

            args = first_quadrant_cloud_extraction.first_quad_dendrogram(min_value=min_val, min_delta=min_delta, min_npix=min_npix)
            d, catalog, header, metadata = args

            for max_tree in max_tree_list:

                for min_tree in min_tree_list:

                    cloud_catalog = first_quadrant_cloud_extraction.export_firstquad_catalog(args, max_descendants=max_tree, min_descendants=min_tree)

                    run_metrics_on(cloud_catalog, d, min_val, min_delta, max_tree, min_tree, min_npix)

