"""
Calculates the velocity split of a given structure.

The velocity split is defined as the minimum v_cen difference between a structure and its two children.

"""

from __future__ import division

import numpy as np


def calculate_velocity_split(d, catalog):

    velocity_split = np.zeros(len(catalog))

    for i in range(len(catalog)):

        idx = catalog['_idx'][i]
        vcen = catalog['v_cen'][i]

        struct = d[idx]

        try:
            child1, child2 = struct.children
        except ValueError:
            # velocity_split[i] gets zero by default, so just move on
            continue

        idx1, idx2 = child1.idx, child2.idx
        vcen1 = catalog['v_cen'][catalog['_idx']==idx1]
        vcen2 = catalog['v_cen'][catalog['_idx']==idx2]

        delta1 = np.abs(vcen - vcen1)
        delta2 = np.abs(vcen - vcen2)

        split = min(delta1, delta2)

        velocity_split[i] = split

    return velocity_split


def descendants_max_vsplit(d, catalog):
    """
    Assigns each structure the biggest v_split of itself or its descendants

    """

    unchecked = -1

    max_vsplit = np.ones(len(catalog)) * unchecked

    # if you're a single feature you get zero
    max_vsplit[catalog['n_descendants'] == 0] = 0

    for i, struct in enumerate(d):

        if max_vsplit[struct.idx] == unchecked:

            # probably hella slow. IDC.
            descendant_idx_list = [x.idx for x in struct.descendants] + [struct.idx]
            descendant_vsplit_list = catalog['v_split'][np.in1d(catalog['_idx'], descendant_idx_list)]

            max_vsplit[struct.idx] = max(descendant_vsplit_list)

    return max_vsplit



