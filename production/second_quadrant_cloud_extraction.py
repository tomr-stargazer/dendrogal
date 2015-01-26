""" Production code to extract a cloud catalog from the second quadrant. """

import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from .load_and_process_data import load_data, permute_data_to_standard_order
from .compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog
from .calculate_distance_dependent_properties import assign_properties
from .remove_degenerate_structures import reduce_catalog

from ..reid_distance_assigner import make_reid_distance_column

def second_quad_dendrogram():
    datacube, header = permute_data_to_standard_order(*load_data("DHT17_Quad2_bw_mominterp.fits"))
    d = compute_dendrogram(datacube, header, min_value=0.31, min_delta=0.31, min_npix=20)
    catalog, metadata = compute_catalog(d, header)

    # DISTANCE assignment
    reid_distance = make_reid_distance_column(catalog, nearfar='near')
    catalog['distance'] = reid_distance['D_k']

    # assignment of physical properties
    assign_properties(catalog)

    return d, catalog, header, metadata


def disqualify(input_catalog):

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (np.abs(catalog['v_cen']) < 25) |
        (catalog['major_sigma'] > 2 * u.deg) |
        (catalog['area_exact'] > 50 * u.deg**2) |
        (catalog['mass'] < 10**3.5 * u.solMass) )

    output_catalog = catalog[~disqualified]

    return output_catalog

def prune_catalog(d, catalog):
    new_catalog = disqualify(catalog)

    smaller_catalog = reduce_catalog(d, new_catalog)

    return smaller_catalog



