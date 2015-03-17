""" Production code to extract a cloud catalog from the second quadrant. """

from __future__ import division

import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from .load_and_process_data import load_data, permute_data_to_standard_order
from .compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog
from .calculate_distance_dependent_properties import assign_properties
from .remove_degenerate_structures import reduce_catalog
from .disqualify_edge_structures import identify_edge_structures

from ..reid_distance_assigner import make_reid_distance_column
from ..catalog_tree_stats import compute_tree_stats

def second_quad_dendrogram():
    datacube, header = permute_data_to_standard_order(*load_data("DHT17_Quad2_bw_mominterp.fits"))
    d = compute_dendrogram(datacube, header, min_value=0.31/2, min_delta=0.31, min_npix=20)
    catalog, metadata = compute_catalog(d, header)

    # DISTANCE assignment
    reid_distance = make_reid_distance_column(catalog, nearfar='near')
    catalog['distance'] = reid_distance['D_k']

    # assignment of physical properties
    assign_properties(catalog)

    # assignment of tree statistic properties
    compute_tree_stats(catalog, d)    

    # note edge structures
    catalog['on_edge'] = identify_edge_structures(d)

    return d, catalog, header, metadata


def extract_clouds(input_catalog):

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (np.abs(catalog['v_cen']) < 25) |
        # (catalog['major_sigma'] > 2 * u.deg) |
        # (catalog['area_exact'] > 50 * u.deg**2) |
        (catalog['fractional_gain'] > 0.81) | 
        (catalog['n_descendants'] > 15) |
        (catalog['on_edge'] == 1) | # removes edge objects
        (catalog['x_sol'] < -7.5) | # screens out kdist errors near l=180
        (catalog['v_cen'] < -100 ) | # no parts of the galaxy are going this fast
        (catalog['mass'] < 10**3.5 * u.solMass) )

    output_catalog = catalog[~disqualified]

    return output_catalog

def prune_catalog(d, catalog):
    new_catalog = extract_clouds(catalog)

    smaller_catalog = reduce_catalog(d, new_catalog)

    return smaller_catalog

def export_secondquad_catalog(args=None):
    """ 
    Uses the above functions to create a "polished" and "final" cloud catalog from this quadrant.

    """

    if args is None:
        d, catalog, header, metadata = second_quad_dendrogram()
    else:
        d, catalog, header, metadata = args

    cloud_catalog = extract_clouds(catalog)

    reduced_cloud_catalog = reduce_catalog(d, cloud_catalog)

    return reduced_cloud_catalog



