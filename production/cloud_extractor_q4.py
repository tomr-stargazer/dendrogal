"""
Extracts clouds from the fourth quadrant using the new cloud-first method.

"""

from __future__ import division

import numpy as np

import astropy
import astropy.units as u

from dendrogal.production.calculate_distance_dependent_properties import assign_properties
from dendrogal.production.remove_degenerate_structures import reduce_catalog
from dendrogal.production.disqualify_edge_structures import identify_edge_structures
from dendrogal.production.distance_disambiguate import distance_disambiguator, assign_distance_columns
from dendrogal.production.velocity_split import calculate_velocity_split, descendants_max_vsplit

from dendrogal.reid_distance_assigner import make_reid_distance_column, distance_assigner_with_plusminus_errors
from dendrogal.catalog_tree_stats import compute_tree_stats

from dendrogal.production.make_fourthquad_stub import d, catalog, header, metadata


def fourth_quad_cloud_catalog():
    """
    This generates a 'processed' dendrogram catalog, made
    from the global, pre-loaded `catalog` from the fourth quadrant.

    It identifies edge structures, computes tree statistics,
    determines the near & far distances (but does not resolve them),
    assigns a distance when there is no KDA, and then assigns physical
    properties to the unambiguously-distanced structures.

    Its output is meant to be the input to `get_positive_velocity_clouds()`
    and `get_negative_velocity_clouds()`.

    """

    catalog_cp = catalog.copy(copy_data=True)

    # assignment of tree statistic properties
    compute_tree_stats(catalog_cp, d)

    # note edge structures
    catalog_cp['on_edge'] = identify_edge_structures(d)

    # calculate velocity split
    catalog_cp['v_split'] = calculate_velocity_split(d, catalog_cp)
    catalog_cp['max_vsplit'] = descendants_max_vsplit(d, catalog_cp)

    near_distance_table = make_reid_distance_column(catalog, nearfar='near')
    far_distance_table = make_reid_distance_column(catalog, nearfar='far')

    # DISTANCE assignment
    distance_assigner_with_plusminus_errors(catalog_cp, near_distance_table, distance_column_name='near_distance')
    distance_assigner_with_plusminus_errors(catalog_cp, far_distance_table, distance_column_name='far_distance')

    # where there's no degeneracy, go ahead and apply the thing
    no_degeneracy = catalog_cp['near_distance'] == catalog_cp['far_distance']
    distance_column = np.zeros_like(near_distance_table['D_k']) * np.nan
    error_distance_column_plus = np.zeros_like(near_distance_table['D_k']) * np.nan
    error_distance_column_minus = np.zeros_like(near_distance_table['D_k']) * np.nan
    KDA_resolution_column = np.array(['-']*len(distance_column))

    distance_column[no_degeneracy] = catalog_cp['near_distance'][no_degeneracy]
    error_distance_column_plus[no_degeneracy] = catalog_cp['error_near_distance_plus'][no_degeneracy]
    error_distance_column_minus[no_degeneracy] = catalog_cp['error_near_distance_minus'][no_degeneracy]
    KDA_resolution_column[no_degeneracy] = 'U'

    assign_distance_columns(catalog_cp, distance_column, KDA_resolution_column, 
                            np.zeros_like(distance_column), np.zeros_like(distance_column),
                            error_distance_column_plus, error_distance_column_minus)

    # assignment of physical properties to unambigously-distanced structures
    # let's think critically about whether this step is needed.
    assign_properties(catalog_cp)

    return catalog_cp


def get_negative_velocity_clouds(input_catalog, max_descendants=10):
    """
    Extracts clouds from the negative-velocity region of the fourth quad.

    This is the 4Q's inner Galaxy.
    Here, the KDA applies to most structures, so we first define
    clouds based on (a) position in velocity space, (b) a cloud is not
    on an edge, (c) line widths between 1-10 km/s, (d) a cloud has less
    than `max_descendants` descendants, (e) the `fractional_gain` is
    below 0.81.

    Structures that meet the above criteria go into a pre-candidate-cloud
    list, and that list is "flattened" or "reduced" to remove degenerate
    substructure.

    Then the KDA is disambiguated for the relevant structures using the
    function `distance_disambiguator` and their physical properties are 
    computed.

    A final list of clouds is generated by taking structures with a mass
    greater than 3 x 10^4 solar masses; this is what's returned.

    """

    catalog = input_catalog.copy(copy_data=True)

    # one. grab the clouds we think are real
    disqualified = (
        (catalog['v_cen'] > -3) |
        (catalog['on_edge'] == 1) |
        (catalog['v_rms'] <= 1) |
        (catalog['v_rms'] > 10) |
        (catalog['max_vsplit'] > 3)
    )

    qualified = (
        (catalog['n_descendants'] < max_descendants) &
        (catalog['fractional_gain'] < 0.9))

    pre_output_catalog = catalog[~disqualified & qualified]

    # now it's just got clouds that COULD be real
    almost_output_catalog = reduce_catalog(d, pre_output_catalog)

    # disambiguate distances here
    assign_distance_columns(almost_output_catalog, *distance_disambiguator(almost_output_catalog, ambiguous_threshold=0.001, scatter_in_log10_R=0.265*2))

    assign_properties(almost_output_catalog)

    # now let's do a thing
    final_qualified = (almost_output_catalog['mass'] > 3e4)
    output_catalog = almost_output_catalog[final_qualified]

    return output_catalog


def get_positive_velocity_clouds(input_catalog, max_descendants=10):
    """
    Extracts clouds from the positive-velocity region of the fourth quad.

    This is the 4Q's Outer Galaxy.
    Here, the KDA does not apply. We fourth define clouds based on 
    (a) position in velocity space, 
    (b) a cloud is not on an edge, 
    (c) a cloud has less than `max_descendants` descendants, 
    (d) the `fractional_gain` is below 0.81.

    Structures that meet the above criteria go into a pre-candidate-cloud
    list, and that list is "flattened" or "reduced" to remove degenerate
    substructure and their physical properties are computed.

    A final list of clouds is generated by taking structures with a mass
    greater than 3 x 10^3 solar masses; this is what's returned.

    Note that the mass criterion is an order of magnitude lower here 
    than in the positive-velocity region; this is mainly because
    (a) less confusion means these lower mass clouds are easy to 
        definitely identify, and
    (b) clouds here are overall lower-mass intrinsically, so we'd have a 
        low count if we didn't include them.

    """

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified_location = (
        ((catalog['v_cen'] < 10) & (catalog['x_cen'] > 340)) |
        ((catalog['v_cen'] < 0) & (catalog['x_cen'] < 340)) |
        (catalog['v_cen'] > 100) |
        (catalog['on_edge'] == 1)
    )

    disqualified_tree = (
        (catalog['n_descendants'] > max_descendants) |
        (catalog['fractional_gain'] > 0.9))

    # this step excludes the weird tail near the galactic center
    disqualified_extreme_positive_velocities_near_GC = (
        (catalog['v_cen'] > 10) & (catalog['x_cen'] > 340))

    pre_output_catalog = catalog[~disqualified_tree &
        ~disqualified_location & ~disqualified_extreme_positive_velocities_near_GC]

    # now it's just got clouds that COULD be real
    almost_output_catalog = reduce_catalog(d, pre_output_catalog)

    # these objects already have mass, size etc computed so that's fine
    final_qualified = (almost_output_catalog['mass'] > 3e3)
    output_catalog = almost_output_catalog[final_qualified]

    return output_catalog


def compile_fourthquad_catalog(input_catalog):

    negative_v_catalog = get_negative_velocity_clouds(input_catalog)
    positive_v_catalog = get_positive_velocity_clouds(input_catalog)

    composite_unreduced_catalog = astropy.table.vstack([negative_v_catalog, positive_v_catalog])

    composite_reduced_catalog = reduce_catalog(d, composite_unreduced_catalog)

    return composite_reduced_catalog


def export_fourthquad_catalog():

    return compile_fourthquad_catalog(fourth_quad_cloud_catalog())

