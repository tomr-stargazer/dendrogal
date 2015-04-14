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

    x = distance * np.sin(0.5*np.pi - brad) *np.cos(lrad) - R_0
    x_solar = distance * np.sin(0.5*np.pi - brad) *np.cos(lrad)
    y = distance * np.sin(0.5*np.pi - brad) *np.sin(lrad)
    z = distance * np.cos(0.5*np.pi - brad)

    catalog['x_sol'] = x_solar.to('kpc')
    catalog['y_sol'] = y.to('kpc')
    catalog['z_sol'] = z.to('kpc')

    # soon, make these cyclindrical! R_gal, phi, z.
    catalog['x_gal'] = x.to('kpc')
    catalog['y_gal'] = y.to('kpc')
    catalog['z_gal'] = z.to('kpc')


