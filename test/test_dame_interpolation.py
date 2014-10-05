""" 
Tests for dame_interpolation.py 

should run with py.test 

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose

from ..dame_interpolation import interpolate_spectrum, interpolate_single, interpolate_double, interpolate_datacube

def test_interpolate_single():

	expected = np.arange(10, dtype='float')
	expected[0] = np.nan

	data = np.copy(expected)
	data[5] = np.nan

	assert_allclose(interpolate_single(data), expected)

def test_interpolate_double():

	expected = np.arange(10, dtype='float')
	expected[0] = np.nan

	data = np.copy(expected)
	data[5] = np.nan
	data[6] = np.nan

	assert_allclose(interpolate_double(data), expected)

def test_interpolate_spectrum():

	expected = np.arange(10, dtype='float')
	expected[0] = np.nan

	data = np.copy(expected)
	data[5] = np.nan
	data[6] = np.nan

	data[8] = np.nan

	assert_allclose(interpolate_spectrum(data), expected)

def test_interpolate_datacube():

	expected = np.arange(60, dtype='float').reshape(3,4,5)
	# should not be interpolated ever
	expected[0,0,0] = np.nan

	data = np.copy(expected)

	# these should be interpolated velocity-wise
	data[1,1,2] = np.nan
	data[1,1,3] = np.nan

	data[0,2,1] = np.nan

	# these should be interpolated spatially
	data[2,2,4] = np.nan

	data[1,3,0] = np.nan

	data[0,2,0] = np.nan

	data[1,2,2] = np.nan # will only work if velocity gets interpolated first	

	assert_allclose(interpolate_datacube(data), expected)



