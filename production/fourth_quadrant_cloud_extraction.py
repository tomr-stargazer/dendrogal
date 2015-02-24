""" Production code to extract a cloud catalog from the first quadrant. """

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
from .detect_disparate_distances import detect_disparate_distances

from ..reid_distance_assigner import make_reid_distance_column
from ..catalog_tree_stats import compute_tree_stats

def fourth_quad_dendrogram():
    datacube, header = permute_data_to_standard_order(*load_data("DHT36_Quad4_mominterp.fits"))
    d = compute_dendrogram(datacube, header, min_value=0.12/2, min_delta=0.12/2, min_npix=20)
    catalog, metadata = compute_catalog(d, header)

    # DISTANCE assignment
    best_distance = distance_disambiguator(catalog)
    catalog['distance'] = best_distance

    # assignment of physical properties
    assign_properties(catalog)

    # assignment of tree statistic properties
    compute_tree_stats(catalog, d)

    # note disparate distances
    catalog['disparate'] = detect_disparate_distances(d, catalog)

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

    quad2_fit_constant = 0.48293812090592952
    quad2_fit_power = 0.56796770148326814

    expected_size =  (1/quad2_fit_constant * catalog['v_rms'].data)**(1/quad2_fit_power) * u.pc

    use_near_distance = (np.abs(near_size - expected_size) <= np.abs(far_size - expected_size))
    use_far_distance = (np.abs(near_size - expected_size) > np.abs(far_size - expected_size))

    best_distance[use_near_distance] = near_distance[use_near_distance]
    best_distance[use_far_distance] = far_distance[use_far_distance]

    return best_distance

def extract_positive_velocity_clouds(input_catalog):

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (catalog['v_cen'] < 5) |
        (catalog['x_cen'] > 340) |
        (catalog['mass'] < 10**3.5 * u.solMass) |
        (catalog['disparate'] == 0))

    output_catalog = catalog[~disqualified]

    return output_catalog


def extract_negative_velocity_clouds(input_catalog):
    """
    This is a way to get things in the negative-velocity part of the map, excluding local stuff.
    """

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (catalog['v_cen'] > -5) |
        (catalog['mass'] < 10**3.5 * u.solMass) |
        (catalog['disparate'] == 0)        
        # (np.abs(catalog['fractional_gain'] - 0.5) > 0.05)
        )

    qualified = (
        (catalog['n_descendants'] < 10) & 
        (catalog['n_descendants'] > 1) &
        (catalog['fractional_gain'] < 0.81))

    output_catalog = catalog[~disqualified & qualified]

    return output_catalog

def export_fourthquad_catalog(args=None):
    """ 
    Uses the above functions to create a "polished" and "final" cloud catalog from this quadrant.

    """

    if args is None:
        d, catalog, header, metadata = fourth_quad_dendrogram()
    else:
        d, catalog, header, metadata = args

    negative_v_catalog = extract_negative_velocity_clouds(catalog)
    positive_v_catalog = extract_positive_velocity_clouds(catalog)

    composite_unreduced_catalog = astropy.table.vstack([negative_v_catalog, positive_v_catalog])

    composite_reduced_catalog = reduce_catalog(d, composite_unreduced_catalog)

    return composite_reduced_catalog


