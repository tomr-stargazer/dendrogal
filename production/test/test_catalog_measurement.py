"""
Tests for catalog_measurement.py

should run with py.test from the root of dendrogal

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy.table import Table

from ..catalog_measurement import size_linewidth_slope

def test_size_linewidth_slope():

	expected_exponent = 0.625
	expected_coefficient = 5

	size_column = np.arange(50)
	linewidth_column = expected_coefficient * size_column**expected_exponent

	test_table = Table()
	test_table['size'] = size_column
	test_table['v_rms'] = linewidth_column

	output = size_linewidth_slope(test_table)

	coefficient = output.beta[0]
	exponent = output.beta[1]

	coefficient_uncertainty = output.sd_beta[0]
	exponent_uncertainty = output.sd_beta[1]

	assert_equal(coefficient, expected_coefficient)
	assert_equal(exponent, expected_exponent)

	assert_allclose(coefficient_uncertainty, 0, atol=1e-7)
	assert_allclose(exponent_uncertainty, 0, atol=1e-7)
	

