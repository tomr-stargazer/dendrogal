"""
Tests things in calculate_distance_dependent_properties.py

should run with py.test

"""

from __future__ import division

from numpy.testing import assert_allclose, assert_equal, assert_almost_equal
import numpy as np

import astropy.table
import astropy.units as u

from ..calculate_distance_dependent_properties import compute_galactic_coordinates, assign_size_with_uncertainties

def test_compute_galactic_coordinates():

    R_0 = 8340 * u.pc
    z_0 = 25 * u.pc

    theta = np.arcsin(z_0/R_0)

    l = 90 * u.deg
    b = 0 * u.deg
    d_sun = R_0

    solar_cartesian, gal_cartesian, gal_cylindrical = compute_galactic_coordinates(l, b, d_sun, R_0=R_0)

    x_sol_expected = 0 * u.pc 
    y_sol_expected = R_0
    z_sol_expected = 0 * u.pc

    assert_allclose(solar_cartesian[0], x_sol_expected, atol=1e-7)
    assert_allclose(solar_cartesian[1], y_sol_expected, atol=1e-7)
    assert_allclose(solar_cartesian[2], z_sol_expected, atol=1e-7)

    x_gal_expected = R_0 * np.cos(theta)
    y_gal_expected = -R_0
    z_gal_expected = R_0 * np.sin(theta)

    assert_allclose(gal_cartesian[0], x_gal_expected, atol=1e-7)
    assert_allclose(gal_cartesian[1], y_gal_expected, atol=1e-7)
    assert_allclose(gal_cartesian[2], z_gal_expected, atol=1e-7)

    R_gal_expected = np.sqrt(x_gal_expected**2 + y_gal_expected**2)
    phi_expected = (-45 * u.deg).to(u.rad)

    assert_allclose(gal_cylindrical[0], R_gal_expected, atol=1e-7)
    assert_allclose(gal_cylindrical[1], phi_expected, rtol=1e-3)


def test_size():

    distance = 1 * u.kpc
    error_distance_plus = 0.3 * u.kpc
    error_distance_minus = 0.2 * u.kpc

    radius = 1 * u.deg

    expected_size = 33.16125578789226 * u.pc
    expected_error_size_plus = 9.948376736367678 * u.pc
    expected_error_size_minus = 6.632251157578453 * u.pc

    catalog = astropy.table.Table()
    catalog['distance'] = u.Quantity([distance])
    catalog['error_distance_plus'] = u.Quantity([error_distance_plus])
    catalog['error_distance_minus'] = u.Quantity([error_distance_minus])

    catalog['radius'] = u.Quantity([radius])

    assign_size_with_uncertainties(catalog)

    assert_allclose(catalog['size'][0], expected_size)
    assert_allclose(catalog['error_size_plus'][0], expected_error_size_plus)
    assert_allclose(catalog['error_size_minus'][0], expected_error_size_minus)

