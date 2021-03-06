""" Production code to extract a cloud catalog from the first quadrant. """

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

def first_quad_dendrogram(min_value=0.18,
                          min_delta=0.18,
                          min_npix=20):

    data_filename = "DHT08_Quad1_mominterp.fits"
    dendrogram_kwargs = {'min_value' : min_value,
                         'min_delta' : min_delta,
                         'min_npix' : min_npix}

    d, catalog, header, metadata = load_permute_dendro_catalog(data_filename, **dendrogram_kwargs)

    # DISTANCE assignment
    best_distance = distance_disambiguator(catalog)
    catalog['distance'] = best_distance

    catalog_clipped = catalog.copy(copy_data=True)

    # assignment of physical properties
    assign_properties(catalog)
    assign_properties(catalog_clipped, flux_column_name='flux_clipped')

    catalog['mass_clipped'] = catalog_clipped['mass']
    catalog['virial_alpha_clipped'] = catalog_clipped['virial_alpha']
    catalog['pressure_clipped'] = catalog_clipped['pressure']

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

    quad2_fit_constant = 0.48293812090592952
    quad2_fit_power = 0.56796770148326814

    expected_size =  (1/quad2_fit_constant * catalog['v_rms'].data)**(1/quad2_fit_power) * u.pc

    use_near_distance = (np.abs(near_size - expected_size) <= np.abs(far_size - expected_size))
    use_far_distance = (np.abs(near_size - expected_size) > np.abs(far_size - expected_size))

    use_near_distance_latitude = (np.abs(catalog['y_cen']) > 1)

    best_distance[use_near_distance] = near_distance[use_near_distance]
    best_distance[use_far_distance] = far_distance[use_far_distance]

    #an override based on latitude
    best_distance[use_near_distance_latitude] = near_distance[use_near_distance_latitude]

    return best_distance


def extract_negative_velocity_clouds(input_catalog, max_descendants=30, min_descendants=2):

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (catalog['v_cen'] > -5) |
        (catalog['mass'] < 10**3.5 * u.solMass) | 
        (catalog['disparate'] == 0) | 
        (catalog['on_edge'] == 1)
        )

    # this step excludes the weird tail near the galactic center
    disqualified_extreme_negative_velocities_near_GC = (
        (catalog['v_cen'] < -10) & (catalog['x_cen'] < 20))

    output_catalog = catalog[~disqualified & ~disqualified_extreme_negative_velocities_near_GC]

    return output_catalog

def extract_positive_velocity_clouds(input_catalog, max_descendants=30, min_descendants=2):
    """
    This is a way to get things in the positive-velocity part of the map, excluding local stuff.
    """

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (catalog['v_cen'] < 20) |
#        (catalog['mass'] < 10**3.5 * u.solMass) |
        (catalog['disparate'] == 0) |
        (catalog['on_edge'] == 1)
        # (np.abs(catalog['fractional_gain'] - 0.5) > 0.05)
        )

    qualified_1 = (
        (catalog['n_descendants'] < max_descendants) & 
        (catalog['n_descendants'] >= min_descendants) &
        (catalog['fractional_gain'] < 0.81) & 
        (catalog['mass'] > 10**3.5 * u.solMass))

    qualified_2 = (
        (catalog['n_descendants'] < max_descendants) & 
        (catalog['fractional_gain'] < 0.81) & 
        (catalog['mass'] > 10**5 * u.solMass))

    output_catalog = catalog[~disqualified & (qualified_1|qualified_2)]

    return output_catalog

def extract_positive_velocity_clouds_control(input_catalog):
    """
    This is a way to get things in the positive-velocity part of the map, excluding local stuff.
    """

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (catalog['v_cen'] < 20) |
        (catalog['disparate'] == 0) |
        (catalog['on_edge'] == 1)
        )

    qualified_1 = (
        (catalog['n_descendants'] < 10) & 
        (catalog['n_descendants'] > 1) &
        (catalog['fractional_gain'] < 0.81) & 
        (catalog['mass'] > 10**3.5 * u.solMass))

    # qualified_2 = (
    #     (catalog['n_descendants'] < 10) & 
    #     (catalog['fractional_gain'] < 0.81) & 
    #     (catalog['mass_clipped'] > 10**5 * u.solMass))

    output_catalog = catalog[~disqualified & #(qualified_1|qualified_2)]
                             qualified_1]

    return output_catalog    

def extract_low_velocity_clouds(input_catalog, max_descendants=30, min_descendants=2):
    """
    Get the near-zero velocity Perseus arm clouds.

    """

    catalog = input_catalog.copy(copy_data=True)

    disqualified = (
        (catalog['v_cen'] < -5) |
        (catalog['v_cen'] > 20) |
        (catalog['n_descendants'] > max_descendants) |
        (catalog['mass'] < 10**3.5 * u.solMass) |
        (np.abs(catalog['y_cen'] > 1)) |
        (catalog['on_edge'] == 1)
        )

    qualified = (
        (catalog['distance'] > 9.8) & 
        (catalog['distance'] < 14))

    output_catalog = catalog[~disqualified & qualified]

    return output_catalog


def export_firstquad_catalog(args=None, **kwargs):
    """ 
    Uses the above functions to create a "polished" and "final" cloud catalog from this quadrant.

    """

    if args is None:
        d, catalog, header, metadata = first_quad_dendrogram()
    else:
        d, catalog, header, metadata = args

    negative_v_catalog = extract_negative_velocity_clouds(catalog, **kwargs)
    positive_v_catalog = extract_positive_velocity_clouds(catalog, **kwargs)
    low_v_catalog = extract_low_velocity_clouds(catalog, **kwargs)

    composite_unreduced_catalog = astropy.table.vstack([negative_v_catalog, positive_v_catalog, low_v_catalog])

    composite_reduced_catalog = reduce_catalog(d, composite_unreduced_catalog)

    return composite_reduced_catalog


