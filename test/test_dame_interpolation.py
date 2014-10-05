""" 
Tests for dame_interpolation.py 

should run with py.test 

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose

from ..dame_interpolation import interpolate_spectrum

def test_interpolate_spectrum():

	expected = np.arange(10)

	data = np.copy(expected)
	data[5] = np.nan

	assert_allclose(data, expected)
