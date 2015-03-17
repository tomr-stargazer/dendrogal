"""
Detects which structures are on the edges of some data.

"""

from __future__ import division

import numpy as np

def on_edge(struct, shape):
    """ Checks whether a given struct is on an edge. """

    x_len, y_len, z_len = shape
    x_array, y_array, z_array = struct.indices()

    if ((0 in x_array) or
        (0 in y_array) or
        (0 in z_array) or
        (x_len-1 in x_array) or
        (y_len-1 in y_array) or
        (z_len-1 in z_array)
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
