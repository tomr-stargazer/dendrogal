"""
This is meant to take in a catalog and generate a "Mark Reid"able table.



Terrific. If it's not hard to generate, what we would like is an ascii table with these columns:

ID  l  b  v  mass

Our program currently expects the ID/name to be 13 characters, e.g., 
G045.09+00.32.  The format of the name does not matter, e.g., it could be 
cloud003, but if you could make it 13 characters long by adding 
underscores or something it would save us a little grief at this end. 
I believe the coordinates can have any width and spacing.

"""

from __future__ import division
import os.path

import numpy as np
import astropy.table

data_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/other_catalogs/")


def id_string_maker(cloud_id):
    """
    Appends the correct number of underscores to the start of an ID.

    This is done to comply with the requirements above.

    """

    string_id = str(cloud_id)

    prefix = "_" * (13 - len(string_id))

    assert len(prefix+string_id) == 13

    return prefix+string_id


def make_reid_table(catalog):

    reid_table = astropy.table.Table()

    ids = range(len(catalog))
    string_ids = [id_string_maker(idx) for idx in ids]

    reid_table['ID'] = string_ids
    reid_table['l'] = catalog['x_cen']
    reid_table['b'] = catalog['y_cen']
    reid_table['v'] = catalog['v_cen']
    reid_table['mass'] = catalog['mass']

    return reid_table


def save_reid_table(reid_table):

    reid_table.write(data_path+'reid_table_v1.txt', format='ascii.basic')
