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

    # Case 2: an edge structure towards the "zero" side
    data = np.zeros((3,3,3))
    data[1,1,1] = 1
    data[1,1,0] = 1
    d = Dendrogram.compute(data, min_value=0)

    on_edge = identify_edge_structures(d)

    expected = np.array([1])

    assert_equal(on_edge, expected)

def test_identify_edge_structures_3():

    # Case 3: an edge structure towards the "max" side
    data = np.zeros((3,3,3))
    data[1,1,1] = 1
    data[1,1,2] = 1
    d = Dendrogram.compute(data, min_value=0)

    on_edge = identify_edge_structures(d)

    expected = np.array([1])

    assert_equal(on_edge, expected)

def test_identify_edge_structures_4():

    # Case 4: something more complex.
    data = np.zeros((4,5,6))

    # a non-edge structure
    data[2,2,1] = 1
    data[2,2,2] = 1

    # an edge structure on zeros
    data[0,0,3] = 1
    data[0,1,3] = 1
    data[1,1,3] = 1

    # an edge structure on the max side
    data[1,1,5] = 1
    data[2,1,5] = 1

    d = Dendrogram.compute(data, min_value=0)

    on_edge = identify_edge_structures(d)

    expected_onedge_len = 2
    expected_offedge_len = 1

    assert_equal(len(d), 3)
    assert_equal(np.sum(on_edge), expected_onedge_len)    
    assert_equal(len(on_edge) - np.sum(on_edge), expected_offedge_len)