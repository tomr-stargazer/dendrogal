""" 
Finds out when structures have incompatible distance assignments internally. 

"""

from __future__ import division

import numpy as np

# so the idea is:
# go through *all* the structures in a catalog/dendrogram
# and find out which ones inside it have distances that disagree.

# this would probbaly take like O(n^2) if I don't optimize
# so, first optimization:
# once a struct has been flagged as "disparate", all of its ancestors are too
# and so you can skip them.

# we can probably draw upon the nifty output of catalog_tree_stats to help us.

def detect_disparate_distances(d, catalog):

    if 'distance' not in catalog.colnames or 'n_descendants' not in catalog.colnames:

        raise ValueError("`catalog` must have `distance` and `n_descendants` columns!")

    # a column: shows "not checked yet", "passes", or "fails".

    passes = 1
    fails = 0
    unchecked = 10

    disparate = np.ones(len(catalog), dtype=np.int) * unchecked

    # if you're a single feature you pass
    disparate[catalog['n_descendants'] == 0] = passes

    # for debugging/analysis - remove when code is complete?
    counter=0

    for i, struct in enumerate(d):

        if disparate[struct.idx] == unchecked:

            # THIS DOESN'T WORK AS DESIRED WE HAVE TO FIX IT 
            # Actually, works as desired in 4th quadrant.
            descendant_idx_list = [x.idx for x in struct.descendants]
            descendant_distance_list = catalog['distance'][np.in1d(catalog['_idx'], descendant_idx_list)]

            if max(descendant_distance_list) >= 2*min(descendant_distance_list):
                fail_struct_and_ancestors(struct, disparate, fails=fails)

                # for debugging/analysis - remove when code is complete?
                if len(descendant_idx_list) < 15:
                    print "Struct {0} has {3} descendants. Distance max/min: {1} and {2}".format(
                        struct.idx, max(descendant_distance_list), 
                        min(descendant_distance_list), len(descendant_idx_list))
                    print " (distance ratio: {0:.3f})".format(
                        max(descendant_distance_list)/min(descendant_distance_list))
                    counter+=1
                    if counter > 10:    
                        return

            else:
                disparate[struct.idx] = passes

    return disparate

    # how does the check go?
    # take a struct
    # investigate its descendants' assigned distances
    # if the min & max of those distances are a factor of 2 or more apart, fail this struct


def fail_struct_and_ancestors(struct, disparate_column, fails=0):
    """
    A recursive function that "fails" the ancestors of a failing struct. 

    """

    disparate_column[struct.idx] = fails
    if struct.level > 0:
        fail_struct_and_ancestors(struct.parent, disparate_column, fails=fails)


