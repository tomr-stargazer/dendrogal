"""
To be run with py.test.

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose, assert_equal

import astropy.table

from ..velocity_split import calculate_velocity_split

def test_calculate_velocity_split():

    mock_dendrogram = {}

    class MockStruct():
        def __init__(self, idx, children=[]):
            self.idx = idx
            self.children = children
            mock_dendrogram[idx] = self

    test_struct1 = MockStruct(1)
    test_struct2 = MockStruct(2)
    test_struct3 = MockStruct(3, children=[test_struct1, test_struct2])
    test_struct4 = MockStruct(4)
    test_struct5 = MockStruct(5, children=[test_struct3, test_struct4])

    test_catalog = astropy.table.Table()
    
    test_idx_column = [1, 2, 3, 4, 5]
    test_vcen_column = [-5, 5, -1, 5, 3]
    test_catalog['_idx'] = test_idx_column
    test_catalog['v_cen'] = test_vcen_column

    velocity_split = calculate_velocity_split(mock_dendrogram, test_catalog)

    expected_split = [0, 0, 4, 0, 2]

    assert_equal(velocity_split, expected_split)




