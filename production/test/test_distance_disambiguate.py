"""
Tests for distance_disambiguate.py

should run with py.test from the root of astrodendro_analysis

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_almost_equal

from ..distance_disambiguate import (p_given_sigmas, p_from_size_linewidth,
                                     p_from_latitude, calculate_p_nearfar,
                                     distance_disambiguator)


def test_p_given_sigmas():
    # Values mined from here:
    # https://en.wikipedia.org/wiki/Standard_deviation#Rules_for_normally_distributed_data

    assert_equal(p_given_sigmas(0), 1)
    assert_almost_equal(p_given_sigmas(1), 0.31731, decimal=4)
    assert_almost_equal(p_given_sigmas(-1), 0.31731, decimal=4)
    assert_almost_equal(p_given_sigmas(2), 0.04550, decimal=4)
    assert_almost_equal(p_given_sigmas(-3), 0.00269, decimal=4)
    assert_almost_equal(p_given_sigmas(4), 0.00006334, decimal=4)


def test_p_from_size_linewidth():

    A = 1
    B = 0.5
    scatter_in_log_R = 0.265

    sigma_v = 5

    expected_log_R = 1/B * np.log10(sigma_v/A)
    expected_R = 10**(expected_log_R)

    assert_equal(
        p_from_size_linewidth(expected_R, sigma_v, A, B, scatter_in_log_R), 1)

    R_plus_1logsigma = 10**(expected_log_R + scatter_in_log_R)

    p_plus_1sigma = p_from_size_linewidth(R_plus_1logsigma, sigma_v, A, B, scatter_in_log_R)

    assert_almost_equal(p_plus_1sigma, p_given_sigmas(1))

    R_minus_2logsigma = 10**(expected_log_R - 2*scatter_in_log_R)

    p_minus_2sigma = p_from_size_linewidth(R_minus_2logsigma, sigma_v, A, B, scatter_in_log_R)

    assert_almost_equal(p_minus_2sigma, p_given_sigmas(-2))

def test_p_from_latitude():

    pass