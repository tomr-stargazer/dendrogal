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


def assign_size_with_uncertainties(catalog):
    """
    R = eta sigma_r d

    error(R)/R = error(d)/d

    """

    eta = 1.9

    sky_radius = eta * u.Quantity(catalog['radius'])
    distance = u.Quantity(catalog['distance'])
    error_distance_plus = u.Quantity(catalog['error_distance_plus'])
    error_distance_minus = u.Quantity(catalog['error_distance_minus'])

    size = sky_radius.to(u.rad).value * distance
    error_size_plus = size * (error_distance_plus/distance)
    error_size_minus = size * (error_distance_minus/distance)

    catalog['size'] = size.to(u.pc)
    catalog['error_size_plus'] = error_size_plus.to(u.pc)
    catalog['error_size_minus'] = error_size_minus.to(u.pc)


def assign_mass_with_uncertainties(catalog, flux_column_name='flux_true'):
    """
    M = 4.4 X_2 F d^2

    error(M)/M = [ (error(X_2)/X)^2 + (2 error(d)/d)^2 ]^(1/2)

    """

    distance = u.Quantity(catalog['distance'])
    error_distance_plus = u.Quantity(catalog['error_distance_plus'])
    error_distance_minus = u.Quantity(catalog['error_distance_minus'])

    flux = u.Quantity(catalog[flux_column_name])

    # mass and alpha calculations from http://adsabs.harvard.edu/abs/2008ApJ...679.1338R

    luminosity = flux * distance**2 / (1 * u.steradian)
    X2 = 1.0
    error_X2 = 0.3 # error on X_CO is around +/- 30% according to Bolatto et al. 2013
    mass = 4.4 * X2 * luminosity.to(u.K * u.km/u.s * u.pc**2).value * u.solMass
    error_mass_plus = mass * ((error_X2/X2)**2 + (2 * error_distance_plus/distance)**2)**(1/2)
    error_mass_minus = mass * ((error_X2/X2)**2 + (2 * error_distance_minus/distance)**2)**(1/2)

    catalog['mass'] = mass.to(u.solMass)
    catalog['error_mass_plus'] = error_mass_plus.to(u.solMass)
    catalog['error_mass_minus'] = error_mass_minus.to(u.solMass)


def assign_alpha_with_uncertainties(catalog, flux_column_name='flux_true'):
    """
    alpha = 5 eta sigma_v^2 R / (G M)

    error(alpha)/alpha = [ (error(X_2)/X)^2 + (error(d)/d)^2 ]^(1/2)

    """

    distance = u.Quantity(catalog['distance'])
    error_distance_plus = u.Quantity(catalog['error_distance_plus'])
    error_distance_minus = u.Quantity(catalog['error_distance_minus'])    

    flux = u.Quantity(catalog[flux_column_name])
    sigma_v = u.Quantity(catalog['v_rms'])

    eta = 1.9

    sky_radius = eta * u.Quantity(catalog['radius'])
    size = sky_radius.to(u.rad).value * distance

    # mass and alpha calculations from http://adsabs.harvard.edu/abs/2008ApJ...679.1338R

    luminosity = flux * distance**2 / (1 * u.steradian)
    X2 = 1.0
    error_X2 = 0.3 # error on X_CO is around +/- 30% according to Bolatto et al. 2013

    mass = 4.4 * X2 * luminosity.to(u.K * u.km/u.s * u.pc**2).value * u.solMass    

    virial_parameter = 5 * eta * sigma_v**2 * size / (mass * c.G)
    error_virial_plus = virial_parameter * ((error_X2/X2)**2 + (error_distance_plus/distance)**2)**(1/2)
    error_virial_minus = virial_parameter * ((error_X2/X2)**2 + (error_distance_minus/distance)**2)**(1/2)

    catalog['virial_alpha'] = virial_parameter.decompose()
    catalog['error_virial_alpha_plus'] = error_virial_plus.decompose()
    catalog['error_virial_alpha_minus'] = error_virial_minus.decompose()


def assign_pressure_with_uncertainties(catalog, flux_column_name='flux_true'):
    """
    P = M sigma_v^2 / (4/3 pi R^3)

    error(P)/P = [ (error(X_2)/X)^2 + (error(d)/d)^2 ]^(1/2)


    """

    distance = u.Quantity(catalog['distance'])
    error_distance_plus = u.Quantity(catalog['error_distance_plus'])
    error_distance_minus = u.Quantity(catalog['error_distance_minus'])    

    flux = u.Quantity(catalog[flux_column_name])
    sigma_v = u.Quantity(catalog['v_rms'])

    eta = 1.9

    sky_radius = eta * u.Quantity(catalog['radius'])
    distance = u.Quantity(catalog['distance'])

    size = sky_radius.to(u.rad).value * distance

    luminosity = flux * distance**2 / (1 * u.steradian)
    X2 = 1.0
    error_X2 = 0.3 # error on X_CO is around +/- 30% according to Bolatto et al. 2013

    mass = 4.4 * X2 * luminosity.to(u.K * u.km/u.s * u.pc**2).value * u.solMass    

    pressure = mass * sigma_v**2 / (4/3 * np.pi * size**3)
    pressure_per_k = pressure / c.k_B

    error_pressure_plus = pressure_per_k * ((error_X2/X2)**2 + (error_distance_plus/distance)**2)**(1/2)
    error_pressure_minus = pressure_per_k * ((error_X2/X2)**2 + (error_distance_minus/distance)**2)**(1/2)

    catalog['pressure'] = pressure_per_k.to(u.K * u.cm**-3)
    catalog['error_pressure_plus'] = error_pressure_plus.to(u.K * u.cm**-3)
    catalog['error_pressure_minus'] = error_pressure_minus.to(u.K * u.cm**-3)


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

    # eta = 1.9

    # sky_radius = eta * u.Quantity(catalog['radius'])
    # distance = u.Quantity(catalog['distance'])

    # size = sky_radius.to(u.rad).value * distance

    # catalog['size'] = size.to(u.pc)
    assign_size_with_uncertainties(catalog)

    # flux = u.Quantity(catalog[flux_column_name])
    # sigma_v = u.Quantity(catalog['v_rms'])

    # # mass and alpha calculations from http://adsabs.harvard.edu/abs/2008ApJ...679.1338R

    # luminosity = flux * distance**2 / (1 * u.steradian)
    # X2 = 1.0
    # mass = 4.4 * X2 * luminosity.to(u.K * u.km/u.s * u.pc**2).value * u.solMass

    # catalog['mass'] = mass.to(u.solMass)
    assign_mass_with_uncertainties(catalog)

    # virial_parameter = 5 * eta * sigma_v**2 * size / (mass * c.G)

    # catalog['virial_alpha'] = virial_parameter.decompose()
    assign_alpha_with_uncertainties(catalog)

    # # Pressure is calculated as a kinetic energy density: mass density times velocity squared
    # pressure = mass * sigma_v**2 / (4/3 * np.pi * size**3)
    # pressure_per_k = pressure / c.k_B

    # catalog['pressure'] = pressure_per_k.to(u.K * u.cm**-3)
    assign_pressure_with_uncertainties(catalog)

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
