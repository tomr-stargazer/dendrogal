"""
Detects which structures are on the edges of some data.

"""

from __future__ import division

import numpy as np

def on_edge(struct, shape):
    """ Checks whether a given struct is on an edge. """

    x_len, y_len, z_len = shape
    x_array, y_array, z_array = struct.indices()

    if ((x_array.min() == 0) or
        (y_array.min() == 0) or
        (z_array.min() == 0) or
        (x_array.max() == x_len-1) or
        (y_array.max() == y_len-1) or
        (z_array.max() == z_len-1)
        ):
        return True
    else:
        return False

def identify_edge_structures(d):
    """ 
    Creates an array identifying which structures are on an edge. 

    Sorted by _idx, not by height - so it follows a catalog order, 
    not a dendrogram order.

    """

    shape = d.data.shape
    on_edge_array = np.zeros(len(d), dtype=int)

    for i in range(len(d)):
        on_edge_array[i] = int(on_edge(d[i], shape))

    return on_edge_array
