"""
Tests functions in remove_degenerate_structures.py

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astrodendro import Structure

from ..remove_degenerate_structures import remove_degenerate_structures

def test_remove_degenerate_structures():

	mock_indices = np.array([0])
	mock_values = np.array([1])

	child3 = Structure(mock_indices, mock_values, idx=5)
	child4 = Structure(mock_indices, mock_values, idx=6)

	child1 = Structure(mock_indices, mock_values, children=[child3, child4], idx=3)
	child2 = Structure(mock_indices, mock_values, idx=4)

	base1 = Structure(mock_indices, mock_values, children=[child1, child2], idx=1)
	base2 = Structure(mock_indices, mock_values, idx=2)

	list_of_structs = [base1, base2, child1, child2, child3, child4]

	expected_list = [base1, base2]

	reduced_list = remove_degenerate_structures(list_of_structs)

	assert_equal(set(reduced_list), set(expected_list))



