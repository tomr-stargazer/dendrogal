""" 
Here is a function that reduces a selection to its principal branches.

"""

import numpy as np

def remove_degenerate_structures(structures):
    """
    Remove all degenerate structures in a list 

    (i.e., structures that are descendants of other structures) 

    """

    structure_set = set(structures)

    # sort by size
    structures_by_descendants = sorted(structures, key=lambda k: len(k.descendants))

    # from smallest to largest:
    for little_struct, i in zip(structures_by_descendants, range(len(structures_by_descendants))):

        # ask from largest to smallest:
        for big_struct in structures_by_descendants[:i:-1]: # backwards, and skip everyone smaller than little_struct

            # "are you my ancestor?" (# wait, why isn't there an "ancestors" property?)
            if little_struct in big_struct.descendants: 

                # if so, remove the little guy
                structure_set.remove(little_struct)
                # and stop asking
                break 

    return list(structure_set)

def reduce_catalog(d, catalog):

    # first, get a list of structures by harvesting a catalog
    struct_list = []
    for row in catalog:
        struct_list.append(d[row['_idx']])

    smaller_struct_list = remove_degenerate_structures(struct_list)

    # now harvest the id's and rebuild the catalog
    smaller_idx_list = [struct.idx for struct in smaller_struct_list]

    smaller_catalog = catalog[np.in1d(catalog['_idx'], smaller_idx_list)]

    return smaller_catalog

def catalog_from_selection(struct_list, catalog):

    selected_ids = [struct.idx for struct in struct_list]
    selected_catalog = catalog[ np.in1d(catalog['_idx'], selected_ids)]

    return selected_catalog

def selection_from_catalog(d, catalog):

    struct_list = [d[idx] for idx in catalog['_idx']]

    return struct_list