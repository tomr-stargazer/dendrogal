""" Production code to extract a cloud catalog from the third quadrant. """

from __future__ import division

import numpy as np

import astropy.units as u

from .convenience_function import load_permute_dendro_catalog
from .calculate_distance_dependent_properties import assign_properties
from .remove_degenerate_structures import reduce_catalog
from .disqualify_edge_structures import identify_edge_structures

from ..reid_distance_assigner import make_reid_distance_column
from ..catalog_tree_stats import compute_tree_stats

def third_quad_dendrogram():

    data_filename = "DHT31_Quad3_mominterp.fits"
    dendrogram_kwargs = {'min_value' : 0.12/2,
                         'min_delta' : 0.12,
                         'min_npix' : 20}

    d, catalog, header, metadata = load_permute_dendro_catalog(data_filename, **dendrogram_kwargs)

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
        (catalog['on_edge'] == 1) | # removes edge objects
        # (catalog['major_sigma'] > 2 * u.deg) |
        # (catalog['area_exact'] > 50 * u.deg**2) |
        (catalog['x_sol'] < -7) |
        (catalog['mass'] < 10**3.5 * u.solMass) )

    output_catalog = catalog[~disqualified]

    return output_catalog

def prune_catalog(d, catalog):
    new_catalog = extract_clouds(catalog)

    smaller_catalog = reduce_catalog(d, new_catalog)

    return smaller_catalog

def export_thirdquad_catalog(args=None):
    """ 
    Uses the above functions to create a "polished" and "final" cloud catalog from this quadrant.

    """

    if args is None:
        d, catalog, header, metadata = third_quad_dendrogram()
    else:
        d, catalog, header, metadata = args

    cloud_catalog = extract_clouds(catalog)

    reduced_cloud_catalog = reduce_catalog(d, cloud_catalog)

    return reduced_cloud_catalog



