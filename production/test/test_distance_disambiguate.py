"""
Tests for distance_disambiguate.py

should run with py.test from the root of dendrogal

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_almost_equal

import astropy.units as u
import astropy.table

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

    # test the most basic functionality
    HWHM = 60 * u.pc

    assert_equal(p_from_latitude(0, HWHM), 1)
    assert_equal(p_from_latitude(-HWHM/np.sqrt(2*np.log(2)), HWHM), p_given_sigmas(1))

    # test, among other things, that for a given latitude and 2 distances,
    # the near distance is more likely

    latitude = 0.5 * u.deg
    near_distance = 5 * u.kpc
    far_distance = 15 * u.kpc
    z_near = near_distance * np.cos(0.5*np.pi - latitude.to(u.rad).value)
    z_far = far_distance * np.cos(0.5*np.pi - latitude.to(u.rad).value)

    p_near = p_from_latitude(z_near, HWHM)
    p_far = p_from_latitude(z_far, HWHM)

    assert p_near > p_far


def test_distance_disambiguator():
    """
    Makes sure that we can determine near/far/ambiguous properly.

    """

    mock_catalog = astropy.table.Table()

    mock_catalog['distance'] = u.Quantity([np.nan]*3, unit=u.kpc)
    mock_catalog['error_distance_plus'] = u.Quantity([np.nan]*3, unit=u.kpc)
    mock_catalog['error_distance_minus'] = u.Quantity([np.nan]*3, unit=u.kpc)

    # the idea is to test that one of them resolves "near" and one "far"
    mock_catalog['near_distance'] = u.Quantity([2]*3, unit=u.kpc)
    mock_catalog['error_near_distance_plus'] = u.Quantity([0.1]*3, unit=u.kpc)
    mock_catalog['error_near_distance_minus'] = u.Quantity([0.1]*3, unit=u.kpc)
    mock_catalog['far_distance'] = u.Quantity([10]*3, unit=u.kpc)
    mock_catalog['error_far_distance_plus'] = u.Quantity([0.1]*3, unit=u.kpc)
    mock_catalog['error_far_distance_minus'] = u.Quantity([0.1]*3, unit=u.kpc)

    # this whole mess is overly coupled...
    mock_catalog['flux_true'] = [1]*3
    mock_catalog['flux_true'].unit = 'K km sr / s'

    # remember R = 1.9 sigma_r d
    mock_catalog['radius'] = [0.015, 0.003, 0.1]
    mock_catalog['radius'].unit = u.deg

    mock_catalog['v_rms'] = [1, 1, 1000]
    mock_catalog['v_rms'].unit = u.km/u.s

    mock_catalog['x_cen'] = [30]*3
    mock_catalog['y_cen'] = [0, 0, 1]

    A_coefficient=1
    B_coefficient=1

    (best_distance, KDA_resolution, p_near, p_far, 
        error_best_distance_plus, 
        error_best_distance_minus) = distance_disambiguator(
        mock_catalog, A_coefficient=A_coefficient, B_coefficient=B_coefficient)

    expected_KDA_resolution = ['N', 'F', 'A']
    expected_distance = [2, 10, np.nan]

    assert_equal(KDA_resolution, expected_KDA_resolution)
    assert_equal(best_distance, expected_distance)