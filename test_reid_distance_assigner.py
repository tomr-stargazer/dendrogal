"""
To be run with py.test.

"""

from __future__ import division

import pdb

import astropy.table
import numpy as np

from .reid_distance_assigner import distance_assigner_with_plusminus_errors


def test_distance_assigner_with_plusminus_errors():

    test_structure_table = astropy.table.Table()

    test_kdist_catalog = astropy.table.Table()
    test_kdist_catalog['D_k'] = np.ones(10)
    test_kdist_catalog['error_D_k_plus'] = np.ones(10)
    test_kdist_catalog['error_D_k_minus'] = -np.ones(10)

    expected_structure_table = astropy.table.Table()
    expected_structure_table['distance'] = np.ones(10)
    expected_structure_table['error_distance_plus'] = np.ones(10)
    expected_structure_table['error_distance_minus'] = -np.ones(10)    

    distance_assigner_with_plusminus_errors(test_structure_table, test_kdist_catalog, distance_column_name='distance')

    pdb.set_trace()

    assert (expected_structure_table == test_structure_table).all()