"""
Tests for disqualify_edge_structures.py

should run with py.test

"""

from numpy.testing import assert_allclose, assert_equal
import numpy as np

from astrodendro import Dendrogram

from ..disqualify_edge_structures import identify_edge_structures

def test_identify_edge_structures_1():

    # Case 1: no edge structures
    data = np.zeros((3,3,3))
    data[1,1,1] = 1
    d = Dendrogram.compute(data, min_value=0)

    on_edge = identify_edge_structures(d)

    expected = np.array([0])

    assert_equal(on_edge, expected)

def test_identify_edge_structures_2():

    # Case 2: an edge structure
    data = np.zeros((3,3,3))
    data[1,1,1] = 1
    data[1,1,0] = 1
    d = Dendrogram.compute(data, min_value=0)

    on_edge = identify_edge_structures(d)

    expected = np.array([1])

    assert_equal(on_edge, expected)

