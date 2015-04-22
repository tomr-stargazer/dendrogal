""" Alternative approach to define clouds solely from the data. """

from __future__ import division

import numpy as np

import astropy
import astropy.units as u

from dendrogal.production.convenience_function import load_permute_dendro_catalog
from dendrogal.production.calculate_distance_dependent_properties import assign_properties
from dendrogal.production.remove_degenerate_structures import reduce_catalog
from dendrogal.production.detect_disparate_distances import detect_disparate_distances
from dendrogal.production.disqualify_edge_structures import identify_edge_structures

from dendrogal.reid_distance_assigner import make_reid_distance_column
from dendrogal.catalog_tree_stats import compute_tree_stats

from make_firstquad_stub import d, catalog, header, metadata


def first_quad_cloud_catalog():
    """
    This generates a 'processed' dendrogram catalog, made
    from the global, pre-loaded `catalog` from the first quadrant.

    It identifies edge structures, computes tree statistics,
    determines the near & far distances (but does not resolve them),
    assigns a distance when there is no KDA, and then assigns physical
    properties to the unambiguously-distanced structures.

    Its output is meant to be the input to `get_positive_velocity_clouds()`,
    `get_perseus_arm_low_velocity_clouds()`, and
    `get_negative_velocity_clouds()`.

    """

    catalog_cp = catalog.copy(copy_data=True)

    # assignment of tree statistic properties
    compute_tree_stats(catalog_cp, d)

    # note edge structures
    catalog_cp['on_edge'] = identify_edge_structures(d)

    near_distance_table = make_reid_distance_column(catalog, nearfar='near')
    far_distance_table = make_reid_distance_column(catalog, nearfar='far')

    near_distance_column = near_distance_table['D_k']
    far_distance_column = far_distance_table['D_k']

    # DISTANCE assignment
    catalog_cp['near_distance'] = near_distance_column
    catalog_cp['far_distance'] = far_distance_column

    # where there's no degeneracy, go ahead and apply the thing
    no_degeneracy = catalog_cp['near_distance'] == catalog_cp['far_distance']
    distance_column = np.zeros_like(near_distance_column) * np.nan
    distance_column[no_degeneracy] = catalog_cp['near_distance'][no_degeneracy]

    catalog_cp['distance'] = distance_column

    # assignment of physical properties to unambigously-distanced structures
    assign_properties(catalog_cp)

    # note disparate distances
    # catalog_cp['disparate'] = detect_disparate_distances(d, catalog_cp)

    return catalog_cp


def get_positive_velocity_clouds(input_catalog, max_descendants=30):
    """
    This extracts clouds from the positive-velocity region of the first quad.

    Here, the KDA applies to most structures, so we first define
    clouds based on (a) position in velocity space, (b) a cloud is not
    on an edge, (c) line widths between 1-10 km/s, (d) a cloud has less
    than `max_descendants` descendants, (e) the `fractional_gain` is 
    below 0.81.

    Structures that meet the above criteria go into a pre-candidate-cloud
    list, and that list is "flattened" or "reduced" to remove degenerate
    substructure.

    Then the KDA is disambiguated for the relevant structures using the
    function `distance_disambiguator` and 

    """


    catalog = input_catalog.copy(copy_data=True)

    # one. grab the clouds we think are real
    disqualified = (
        (catalog['v_cen'] < 20) |
        #        (catalog['mass'] < 10**3.5 * u.solMass) |
        # (catalog['disparate'] == 0) |
        (catalog['on_edge'] == 1) |
        (catalog['v_rms'] <= 1) |
        (catalog['v_rms'] > 10)
        # (np.abs(catalog['fractional_gain'] - 0.5) > 0.05)
    )

    qualified = (
        (catalog['n_descendants'] < max_descendants) &
        (catalog['fractional_gain'] < 0.81))  # &
    # (catalog['mass'] > 10**5 * u.solMass))

    pre_output_catalog = catalog[~disqualified & qualified]

    # now it's just got clouds that COULD be real
    almost_output_catalog = reduce_catalog(d, pre_output_catalog)

    # disambiguate distances here
    best_distance = distance_disambiguator(almost_output_catalog)
    almost_output_catalog['distance'] = best_distance
    assign_properties(almost_output_catalog)

    # now let's do a thing
    final_qualified = (almost_output_catalog['mass'] > 3e4)
    output_catalog = almost_output_catalog[final_qualified]

    return output_catalog


def distance_disambiguator(catalog):
    near_distance_column = catalog['near_distance']
    far_distance_column = catalog['far_distance']

    best_distance = np.zeros_like(near_distance_column)

    sky_radius = u.Quantity(catalog['radius'].data * catalog['radius'].unit)
    near_distance = u.Quantity(near_distance_column)
    near_size = sky_radius.to(u.rad).value * near_distance

    far_distance = u.Quantity(far_distance_column)
    far_size = sky_radius.to(u.rad).value * far_distance

    quad2_fit_constant = 0.48293812090592952
    quad2_fit_power = 0.56796770148326814

    expected_size = (
        1 / quad2_fit_constant * catalog['v_rms'].data) ** (1 / quad2_fit_power) * u.pc

    use_near_distance = (np.abs(near_size - expected_size) <= np.abs(far_size - expected_size))
    use_far_distance = (np.abs(near_size - expected_size) > np.abs(far_size - expected_size))

    use_near_distance_latitude = (np.abs(catalog['y_cen']) > 1)

    best_distance[use_near_distance] = near_distance[use_near_distance]
    best_distance[use_far_distance] = far_distance[use_far_distance]

    # an override based on latitude
    best_distance[use_near_distance_latitude] = near_distance[use_near_distance_latitude]

    return best_distance


def extract_negative_velocity_clouds(input_catalog, max_descendants=30, min_descendants=2):

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified = (
        (catalog['v_cen'] > -5) |
        (catalog['mass'] < 10 ** 3.5 * u.solMass) |
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
        (catalog['mass'] > 10 ** 3.5 * u.solMass))

    qualified_2 = (
        (catalog['n_descendants'] < max_descendants) &
        (catalog['fractional_gain'] < 0.81) &
        (catalog['mass'] > 10 ** 5 * u.solMass))

    output_catalog = catalog[~disqualified & (qualified_1 | qualified_2)]

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
        (catalog['mass'] > 10 ** 3.5 * u.solMass))

    # qualified_2 = (
    #     (catalog['n_descendants'] < 10) &
    #     (catalog['fractional_gain'] < 0.81) &
    #     (catalog['mass_clipped'] > 10**5 * u.solMass))

    output_catalog = catalog[~disqualified &  # (qualified_1|qualified_2)]
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
        (catalog['mass'] < 10 ** 3.5 * u.solMass) |
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

    composite_unreduced_catalog = astropy.table.vstack(
        [negative_v_catalog, positive_v_catalog, low_v_catalog])

    composite_reduced_catalog = reduce_catalog(d, composite_unreduced_catalog)

    return composite_reduced_catalog
