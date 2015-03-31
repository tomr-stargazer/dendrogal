"""
Uses a grid of cloud extraction parameters to optimize the first quadrant cloud extraction parameters.

"""

from __future__ import division

import numpy as np

from astrodendro_analysis.production import first_quadrant_cloud_extraction, map_figures

max_tree_list = [10, 20, 30]
min_val_list = [0.18/2, 0.18]
min_delta_list = [0.18/2, 0.18]
min_tree_list = [1, 0]

print "Total upcoming iterations: {0}".format(len(max_tree_list)*len(min_val_list)*len(min_delta_list)*len(min_tree_list))

directory = 'grid_of_cloud_parameters/'

def run_metrics_on(cloud_catalog, dendrogram, 
                   min_val=None, min_delta=None, 
                   max_tree=None, min_tree=None):

    # make and save some plots

    fig = map_figures.make_lv_map_new(cloudy1, d1, ellipse_color=colorbrewer_blue)
    min_value_str = '{:.3f}'.format(min_value)
    min_delta_str = '{:.3f}'.format(min_delta)
    max_tree_str = str(max_tree)
    min_tree_str = str(min_tree)
    
    filename = "{0}_{1}_{2}_{3}.png".format(min_value_str, min_delta_str, max_tree_str, min_tree_str)

    fig.savefig(directory+filename, bbox_inches='tight')

for min_val in min_val_list:

    for min_delta in min_delta_list:

        args = first_quadrant_cloud_extraction.first_quad_dendrogram(min_value=min_val, min_delta=min_delta)
        d, catalog, header, metadata = args

        for max_tree in max_tree_list:

            for min_tree in min_tree_list:

                cloud_catalog = export_firstquad_catalog(args, max_tree=max_tree, min_tree=min_tree)

                run_metrics_on(cloud_catalog, d, min_val, min_delta, max_tree, min_tree)

