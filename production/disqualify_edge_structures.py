"""
Detects which structures are on the edges of some data.

"""

from __future__ import division

import numpy as np

def identify_edge_structures(d):

    x_len, y_len, z_len = d.data.shape

    on_edge = np.zeros(len(d), dtype=int)

    for i, struct in enumerate(d):

        x_array, y_array, z_array = struct.indices()
        if ((0 in x_array) or
            (0 in y_array) or
            (0 in z_array) or
            (x_len-1 in x_array) or
            (y_len-1 in y_array) or
            (z_len-1 in z_array)
            ):
            on_edge[i] = 1
        else:
            pass

    return on_edge
