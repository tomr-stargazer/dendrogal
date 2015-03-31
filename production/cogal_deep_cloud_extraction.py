""" Production code to extract clouds in the gaps between quadrants """

from __future__ import division

import numpy as np

import astropy
import astropy.units as u

from .convenience_function import load_permute_dendro_catalog
from .calculate_distance_dependent_properties import assign_properties
from .remove_degenerate_structures import reduce_catalog
from .detect_disparate_distances import detect_disparate_distances
from .disqualify_edge_structures import identify_edge_structures

from ..reid_distance_assigner import make_reid_distance_column
from ..catalog_tree_stats import compute_tree_stats

def cogal_deep_dendrogram():

    print "Warning: the main issue with this code is that the WCS object needs special treatment in order to work properly."

    data_filename = "COGAL_deep_mominterp.fits"
    dendrogram_kwargs = {'min_value' : 0.5,
                         'min_delta' : 0.5,
                         'min_npix' : 2000}

    d, catalog, header, metadata = load_permute_dendro_catalog(data_filename, **dendrogram_kwargs)

    # DISTANCE assignment
    best_distance = distance_disambiguator(catalog)
    catalog['distance'] = best_distance

    # assignment of physical properties
    assign_properties(catalog)

    # assignment of tree statistic properties
    compute_tree_stats(catalog, d)

    # note disparate distances
    catalog['disparate'] = detect_disparate_distances(d, catalog)

    # note edge structures
    catalog['on_edge'] = identify_edge_structures(d)

    return d, catalog, header, metadata

def distance_disambiguator(catalog):
    near_distance_table = make_reid_distance_column(catalog, nearfar='near')
    far_distance_table = make_reid_distance_column(catalog, nearfar='far')

    near_distance_column = near_distance_table['D_k']
    far_distance_column = far_distance_table['D_k']

    best_distance = np.zeros_like(near_distance_table['D_k'])

    sky_radius = u.Quantity(catalog['radius'].data * catalog['radius'].unit)
    near_distance = u.Quantity(near_distance_column)
    near_size = sky_radius.to(u.rad).value * near_distance

    far_distance = u.Quantity(far_distance_column)
    far_size = sky_radius.to(u.rad).value * far_distance

    quad2_fit_constant = 0.39
    quad2_fit_power = 0.62

    expected_size =  (1/quad2_fit_constant * catalog['v_rms'].data)**(1/quad2_fit_power) * u.pc

    use_near_distance = (np.abs(near_size - expected_size) <= np.abs(far_size - expected_size))
    use_far_distance = (np.abs(near_size - expected_size) > np.abs(far_size - expected_size))

    best_distance[use_near_distance] = near_distance[use_near_distance]
    best_distance[use_far_distance] = far_distance[use_far_distance]

    return best_distance
