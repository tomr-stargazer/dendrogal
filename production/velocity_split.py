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


