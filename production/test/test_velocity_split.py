"""
To be run with py.test.

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose, assert_equal

import astropy.table

from ..velocity_split import calculate_velocity_split, descendants_max_vsplit

def test_calculate_velocity_split():

    class MockDendrogram():
        def __init__(self):
            self.dict = {}

        def __getitem__(self, key):
            return self.dict[key]

        def __iter__(self):
            return iter(self.dict.values())

    mock_dendrogram = MockDendrogram()

    class MockStruct():
        def __init__(self, idx, children=[]):
            self.idx = idx
            self.children = children
            self.descendants = children # this is meant to be manually overridden in the one necessary case
            mock_dendrogram.dict[idx] = self

    test_struct1 = MockStruct(0)
    test_struct2 = MockStruct(1)
    test_struct3 = MockStruct(2, children=[test_struct1, test_struct2])
    test_struct4 = MockStruct(3)
    test_struct5 = MockStruct(4, children=[test_struct3, test_struct4])
    test_struct5.descendants = [test_struct1, test_struct2, test_struct3, test_struct4]

    test_catalog = astropy.table.Table()
    
    test_idx_column = [0, 1, 2, 3, 4]
    test_vcen_column = [-5, 5, -1, 5, 3]
    test_catalog['_idx'] = test_idx_column
    test_catalog['v_cen'] = test_vcen_column

    velocity_split = calculate_velocity_split(mock_dendrogram, test_catalog)

    expected_split = [0, 0, 4, 0, 2]

    assert_equal(velocity_split, expected_split)

    test_catalog['v_split'] = velocity_split
    test_catalog['n_descendants'] = [0, 0, 2, 0, 4]
    max_vsplit = descendants_max_vsplit(mock_dendrogram, test_catalog)

    excpected_max_split = [0, 0, 4, 0, 4]
    assert_equal(max_vsplit, excpected_max_split)

