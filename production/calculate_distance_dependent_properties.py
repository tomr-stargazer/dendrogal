"""
Takes a catalog and calculates distance-dependent properties.

Most of these relationships are taken from Rosolowsky et al. (2008),
which in turn largely borrowed from Rosolowsky & Leroy (2006).

"""

from __future__ import division

import numpy as np

import astropy
import astrodendro

import astropy.units as u
import astropy.constants as c


def assign_properties(catalog, galactic_center_distance=8.340*u.kpc, flux_column_name='flux_true'):
    """
    Computes and assigns distance-dependent properties to a catalog.

    These properties are:
    size (pc)
    mass (solMass)
    virial_alpha 
    pressure (K cm-3)
    Solar coordinates x, y, z (kpc)
    Galactic coordinates x, y, z (kpc)

    """

    if 'distance' not in catalog.colnames:
        raise ValueError("`catalog` must have a `distance` column")

    eta = 1.9

    sky_radius = eta * u.Quantity(catalog['radius'])
    distance = u.Quantity(catalog['distance'])

    size = sky_radius.to(u.rad).value * distance

    catalog['size'] = size.to(u.pc)

    flux = u.Quantity(catalog[flux_column_name])
    sigma_v = u.Quantity(catalog['v_rms'])

    # mass and alpha calculations from http://adsabs.harvard.edu/abs/2008ApJ...679.1338R

    luminosity = flux * distance**2 / (1 * u.steradian)
    X2 = 1.0
    mass = 4.4 * X2 * luminosity.to(u.K * u.km/u.s * u.pc**2).value * u.solMass

    catalog['mass'] = mass.to(u.solMass)

    virial_parameter = 5 * eta * sigma_v**2 * size / (mass * c.G)

    catalog['virial_alpha'] = virial_parameter.decompose()

    # Pressure is calculated as a kinetic energy density: mass density times velocity squared
    pressure = mass * sigma_v**2 / (4/3 * np.pi * size**3)
    pressure_per_k = pressure / c.k_B

    catalog['pressure'] = pressure_per_k.to(u.K * u.cm**-3)

    R_0 = u.Quantity(galactic_center_distance, u.kpc)

    distance = catalog['distance'].data * catalog['distance'].unit
    lrad = (catalog['x_cen'] * u.deg).to(u.rad).value
    brad = (catalog['y_cen'] * u.deg).to(u.rad).value

    solar_cart, gal_cart, gal_cyl = compute_galactic_coordinates(lrad, brad, distance, R_0=R_0)

    catalog['x_sol'] = solar_cart[0].to(u.kpc)
    catalog['y_sol'] = solar_cart[1].to(u.kpc)
    catalog['z_sol'] = solar_cart[2].to(u.kpc)

    catalog['x_gal'] = gal_cart[0].to(u.kpc)
    catalog['y_gal'] = gal_cart[1].to(u.kpc)
    catalog['z_gal'] = gal_cart[2].to(u.kpc)

    catalog['R_gal'] = gal_cyl[0].to(u.kpc)
    catalog['phi_gal'] = gal_cyl[1].to(u.deg)



def compute_galactic_coordinates(l, b, d_sun, R_0=8.340*u.kpc):
    """
    Cartesian & cylindrical Solar & Galactic coordinates, taking solar offset into account.

    Calculation taken from Ellsworth-Bowers et al. (2013), ApJ 770:39
    Appendix C: "The vertical solar offset and converting (l, b, d_sun)
    into (R_gal, phi, z)"

    """

    from numpy import cos, sin, arcsin

    x_sol = d_sun * cos(l) * cos(b)
    y_sol = d_sun * sin(l) * cos(b)
    z_sol = d_sun * sin(b)

    solar_cartesian = (x_sol, y_sol, z_sol)

    # Sun's displacement above the midplane
    z_0 = 25 * u.pc
    theta = arcsin(z_0/R_0)

    x_gal = R_0*cos(theta) - d_sun * (cos(l)*cos(b)*cos(theta) + sin(b)*sin(theta))
    y_gal = -d_sun * sin(l) * cos(b)
    z_gal = R_0*sin(theta) - d_sun * (cos(l)*cos(b)*sin(theta) + sin(b)*cos(theta))

    gal_cartesian = (x_gal, y_gal, z_gal)

    R_gal = np.sqrt(x_gal**2 + y_gal**2)
    phi_gal = np.arctan2(y_gal, x_gal)
    z_gal = z_gal

    gal_cylindrical = (R_gal, phi_gal, z_gal)

    return solar_cartesian, gal_cartesian, gal_cylindrical
