""" 
Extracts clouds from the Carina survey using the new cloud-first method.

"""

from __future__ import division

import numpy as np

import astropy
import astropy.units as u

from dendrogal.production.convenience_function import load_permute_dendro_catalog
from dendrogal.production.calculate_distance_dependent_properties import assign_properties
from dendrogal.production.remove_degenerate_structures import reduce_catalog
from dendrogal.production.disqualify_edge_structures import identify_edge_structures
from dendrogal.production.distance_disambiguate import distance_disambiguator


from dendrogal.reid_distance_assigner import make_reid_distance_column
from dendrogal.catalog_tree_stats import compute_tree_stats

from make_carina_stub import d, catalog, header, metadata


def carina_cloud_catalog():
    """
    This generates a 'processed' dendrogram catalog, made
    from the global, pre-loaded `catalog` from the Carina survey.

    It identifies edge structures, computes tree statistics,
    determines the near & far distances (but does not resolve them),
    assigns a distance when there is no KDA, and then assigns physical
    properties to the unambiguously-distanced structures.

    Its output is meant to be the input to `get_positive_velocity_clouds()`,
    and `get_negative_velocity_clouds()`.

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
    # let's think critically about whether this step is needed.
    assign_properties(catalog_cp)

    return catalog_cp

def get_negative_velocity_clouds(input_catalog, max_descendants=10):
    """
    Extracts clouds from the negative-velocity region of the Carina survey.

    This is the 4Q's inner Galaxy.
    Here, the KDA might apply to some structures, so we first define
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
        (catalog['v_cen'] > -15) |
        (catalog['on_edge'] == 1) |
        (catalog['v_rms'] <= 1) |
        (catalog['v_rms'] > 10) |
        (catalog['x_cen'] < 280)
    )

    qualified = (
        (catalog['n_descendants'] < max_descendants) &
        (catalog['fractional_gain'] < 0.81))

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

def get_positive_velocity_clouds(input_catalog, max_descendants=30):
    """
    Extracts clouds from the positive-velocity region of Carina.

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

    disqualified_location = (
        (catalog['v_cen'] < 10) |
        (catalog['on_edge'] == 1) |
        (catalog['x_cen'] < 280)        
    )

    disqualified_tree = (
        (catalog['n_descendants'] > max_descendants) &
        (catalog['fractional_gain'] > 0.81))

    pre_output_catalog = catalog[~disqualified_location & ~disqualified_tree]

    # now it's just got clouds that COULD be real
    almost_output_catalog = reduce_catalog(d, pre_output_catalog)

    # these objects already have mass, size etc computed so that's fine
    final_qualified = (almost_output_catalog['mass'] > 3e3)
    output_catalog = almost_output_catalog[final_qualified]

    return output_catalog


def compile_carina_catalog(input_catalog):

    negative_v_catalog = get_negative_velocity_clouds(input_catalog)
    positive_v_catalog = get_positive_velocity_clouds(input_catalog)

    composite_unreduced_catalog = astropy.table.vstack([negative_v_catalog, positive_v_catalog])

    composite_reduced_catalog = reduce_catalog(d, composite_unreduced_catalog)

    return composite_reduced_catalog


def export_carina_catalog():

    return compile_carina_catalog(carina_cloud_catalog())
