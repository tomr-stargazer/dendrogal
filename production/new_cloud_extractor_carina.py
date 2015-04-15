""" Production code to extract a cloud catalog from the Carina survey. """

from __future__ import division

import numpy as np

import astropy.units as u

from astrodendro_analysis.production.convenience_function import load_permute_dendro_catalog
from astrodendro_analysis.production.calculate_distance_dependent_properties import assign_properties
from astrodendro_analysis.production.remove_degenerate_structures import reduce_catalog
from astrodendro_analysis.production.disqualify_edge_structures import identify_edge_structures

from astrodendro_analysis.reid_distance_assigner import make_reid_distance_column
from astrodendro_analysis.catalog_tree_stats import compute_tree_stats


def carina_dendrogram():

    data_filename = "DHT33_Carina_mominterp.fits"
    dendrogram_kwargs = {'min_value' : 0.17/2,
                         'min_delta' : 0.17,
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
        (catalog['mass'] < 10**3.5 * u.solMass) )

    output_catalog = catalog[~disqualified]

    return output_catalog


def export_carina_catalog(args=None):
    """ 
    Uses the above functions to create a "polished" and "final" cloud catalog from this quadrant.

    """

    if args is None:
        d, catalog, header, metadata = carina_dendrogram()
    else:
        d, catalog, header, metadata = args

    cloud_catalog = extract_clouds(catalog)

    reduced_cloud_catalog = reduce_catalog(d, cloud_catalog)

    return reduced_cloud_catalog



