""" A module to compute tree statistics for a catalog. """

from __future__ import division

import numpy as np

def compute_tree_stats(catalog, dendrogram):

    if len(catalog) > len(dendrogram):
        raise ValueError("dendrogram cannot have fewer entries than catalog")

    level_list = []
    n_descendants_list = []
    is_leaf_list = []
    fractional_gain_list = []

    for idx in catalog['_idx']:

        struct = dendrogram[idx]

        level = struct.level
        n_descendants = len(struct.descendants)
        is_leaf = struct.is_leaf
        try:
            fractional_gain = n_descendants / len(struct.parent.descendants)
        except AttributeError:
            fractional_gain = np.nan

        level_list.append(level)
        n_descendants_list.append(n_descendants)
        is_leaf_list.append(is_leaf)
        fractional_gain_list.append(fractional_gain)

    catalog['level'] = level_list
    catalog['n_descendants'] = n_descendants_list
    catalog['is_leaf'] = is_leaf_list
    catalog['fractional_gain'] = fractional_gain_list        

    return None