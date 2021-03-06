"""
Extracts clouds from the second quadrant using the new cloud-first method.

"""

from __future__ import division

import numpy as np

import astropy

from dendrogal.production.calculate_distance_dependent_properties import assign_properties
from dendrogal.production.remove_degenerate_structures import reduce_catalog
from dendrogal.production.disqualify_edge_structures import identify_edge_structures
from dendrogal.production.convenience_function import load_permute_dendro_catalog
from dendrogal.production.distance_disambiguate import assign_distance_columns_trivial

from dendrogal.reid_distance_assigner import make_universal_distance_column, distance_assigner_with_plusminus_errors
from dendrogal.catalog_tree_stats import compute_tree_stats

from dendrogal.production.make_secondquad_stub import d, catalog, header, metadata

def second_quad_cloud_catalog():
    """
    This generates a 'processed' dendrogram catalog, made
    from the global, pre-loaded `catalog` from the second quadrant survey.

    It identifies edge structures, computes tree statistics,
    determines the kinematic distances, and then assigns physical
    properties to the structures.

    Its output is meant to be the input to `get_positive_velocity_clouds()`,
    and `get_negative_velocity_clouds()`.

    """

    catalog_cp = catalog.copy(copy_data=True)

    # assignment of tree statistic properties
    compute_tree_stats(catalog_cp, d)

    # note edge structures
    catalog_cp['on_edge'] = identify_edge_structures(d)

    # DISTANCE assignment
    distance_table = make_universal_distance_column(catalog, nearfar='near')

    distance_assigner_with_plusminus_errors(catalog_cp, distance_table, distance_column_name='distance')
    assign_distance_columns_trivial(catalog_cp)

    # assignment of physical properties to unambigously-distanced structures
    assign_properties(catalog_cp)

    return catalog_cp


def get_negative_velocity_clouds(input_catalog, max_descendants=10):

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified_location = (
        (catalog['v_cen'] > -8) |
        (catalog['on_edge'] == 1) |
        (catalog['v_cen'] < -100) |
        (catalog['x_sol'] < -7.5)
    )

    disqualified_tree = (
        (catalog['n_descendants'] > max_descendants) |
        (catalog['fractional_gain'] > 0.9))


    pre_output_catalog = catalog[~disqualified_tree & ~disqualified_location]

    # now it's just got clouds that COULD be real
    almost_output_catalog = reduce_catalog(d, pre_output_catalog)

    # these objects already have mass, size etc computed so that's fine
    final_qualified = (almost_output_catalog['mass'] > 3e3)
    output_catalog = almost_output_catalog[final_qualified]

    return output_catalog

def get_positive_velocity_clouds(input_catalog, max_descendants=10):

    catalog = input_catalog.copy(copy_data=True)

    # narrow down how we select clouds
    disqualified_location = (
        (catalog['v_cen'] < 12) |
        (catalog['on_edge'] == 1) |
        (catalog['x_sol'] < -7.5) |
        (catalog['x_cen'] < 190)
    )

    disqualified_tree = (
        (catalog['n_descendants'] > max_descendants) |
        (catalog['fractional_gain'] > 0.9))

    pre_output_catalog = catalog[~disqualified_tree & ~disqualified_location]

    # now it's just got clouds that COULD be real
    almost_output_catalog = reduce_catalog(d, pre_output_catalog)

    # these objects already have mass, size etc computed so that's fine
    final_qualified = (almost_output_catalog['mass'] > 3e3)
    output_catalog = almost_output_catalog[final_qualified]

    return output_catalog

def compile_secondquad_catalog(input_catalog):

    negative_v_catalog = get_negative_velocity_clouds(input_catalog)
    positive_v_catalog = get_positive_velocity_clouds(input_catalog)

    composite_unreduced_catalog = astropy.table.vstack([negative_v_catalog, positive_v_catalog])

    composite_reduced_catalog = reduce_catalog(d, composite_unreduced_catalog)

    return composite_reduced_catalog

def export_secondquad_catalog():

    return compile_secondquad_catalog(second_quad_cloud_catalog())

