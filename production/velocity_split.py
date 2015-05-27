"""
Calculates the velocity split of a given structure.

The velocity split is defined as the v_cen difference between a structure's two children.

"""

from __future__ import division

import numpy as np


def calculate_velocity_split(d, catalog):

    velocity_split = np.zeros(len(catalog))

    for i in range(len(catalog)):

        idx = catalog['_idx'][i]

        struct = d[idx]

        try:
            child1, child2 = struct.children
        except ValueError:
            # velocity_split[i] gets zero by default, so just move on
            continue

        idx1, idx2 = child1.idx, child2.idx
        vcen1 = catalog['v_cen'][catalog['_idx']==idx1]
        vcen2 = catalog['v_cen'][catalog['_idx']==idx2]

        split = np.abs(vcen2 - vcen1)

        velocity_split[i] = split

    return velocity_split


