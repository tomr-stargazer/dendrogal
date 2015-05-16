"""
To be run with py.test.

"""

from __future__ import division

import pytest

import astropy.table
import numpy as np

def reid_import_fails():
    try:
        from .reid_distance_assigner import distance_assigner_with_plusminus_errors
        return False
    except ImportError:
        return True

@pytest.mark.skipif(reid_import_fails(),
                    reason="requires Reid kdist tool installed")
def test_distance_assigner_with_plusminus_errors():

    from .reid_distance_assigner import distance_assigner_with_plusminus_errors

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

    assert (expected_structure_table == test_structure_table).all()

    test_structure_table2 = astropy.table.Table()

    expected_structure_table2 = astropy.table.Table()
    expected_structure_table2['near_distance'] = np.ones(10)
    expected_structure_table2['error_near_distance_plus'] = np.ones(10)
    expected_structure_table2['error_near_distance_minus'] = -np.ones(10)    

    distance_assigner_with_plusminus_errors(test_structure_table2, test_kdist_catalog, distance_column_name='near_distance')

    assert (expected_structure_table2 == test_structure_table2).all()

