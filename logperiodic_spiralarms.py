""" 
logperiodic_spiralarms.py

"""

from __future__ import division
from functools import partial

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u


def logperiodic(reference_radius, reference_azimuth, pitch_angle, azimuth):
    """
    Computes a log-periodic spiral curve in polar coordinates.

    Parameters
    ----------
    reference_radius : float or Quantity
        Galactocentric radius `R_ref` at reference angle `beta_ref`.
        default unit: kpc
    reference_azimuth : float or Quantity
        Galactocentric reference azimuth `beta_ref` for `R_ref`.
        default unit: deg
    pitch_angle : float or Quantity
        Spiral pitch angle `psi`.
        default unit: deg
    azimuth : np.ndarray or Quantity
        Galactocentric azimuth `beta`.
        default unit: deg

    Returns
    -------
    spiral_radius : Quantity
        Galactocentric radius `R` of spiral at `beta`.
        units: kpc

    """

    Rref = u.Quantity(reference_radius, u.kpc)
    betaref = u.Quantity(reference_azimuth, u.deg)
    psi = u.Quantity(pitch_angle, u.deg)
    beta = u.Quantity(azimuth, u.deg)

    lnR = np.log(Rref/u.kpc) - (beta - betaref).to(u.rad).value * np.tan(psi)
    spiral_radius = np.exp(lnR) * u.kpc

    return spiral_radius

scutum_arm = partial(logperiodic, 5.0, 27.6, 19.8)
sagittarius_arm = partial(logperiodic, 6.6, 25.6, 6.9)
local_arm = partial(logperiodic, 8.4, 8.9, 12.8)
perseus_arm = partial(logperiodic, 9.9, 14.2, 9.4)
outer_arm = partial(logperiodic, 13, 18.6, 13.8)

scutum_angles = np.arange(3, 101, 0.125)
sagittarius_angles = np.arange(-2, 68, 0.125)
local_angles = np.arange(-8, 27, 0.125)
perseus_angles = np.arange(-21, 88, 0.125)
outer_angles = np.arange(-6, 56, 0.125)

scutum_radii = scutum_arm(scutum_angles)
sagittarius_radii = sagittarius_arm(sagittarius_angles)
local_radii = local_arm(local_angles)
perseus_radii =  perseus_arm( perseus_angles)
outer_radii = outer_arm(outer_angles)

def convert_galactic_polar_to_solar_polar(beta, radii, solar_radius=8.34):

    rsun = u.Quantity(solar_radius, u.kpc)
    radii = u.Quantity(radii, u.kpc)
    beta = u.Quantity(beta, u.deg)

    galactic_x = radii * np.sin(np.radians(beta))
    galactic_y = radii * np.cos(np.radians(beta))

    solar_x = galactic_x
    solar_y = galactic_y - rsun

    solar_d = (solar_x**2 + solar_y**2)**(1/2)
    solar_l = np.arctan2( solar_y, solar_x ).to(u.deg) - 90 * u.deg

    # we constrain solar_d to be positive, so sometimes the angle needs to wrap around halfway
    # solar_l[solar_x >= 0] += 180 * u.deg

    # make negative angles wrap until they're positive
    while (solar_l < 0).any():
        solar_l[solar_l < 0] += 360 * u.deg

    return solar_d, solar_l


def arm_distance(solar_d_array, solar_l_arry, longitude):

    longitude = u.Quantity(longitude, u.deg)

    angular_separation = np.abs(longitude - solar_l_arry)

    if angular_separation.min() > 1*u.deg:
        return np.nan

    return solar_d_array[ angular_separation == angular_separation.min() ][0]



# fig = plt.figure()

# plt.polar(np.radians(scutum_angles), scutum_radii)
# plt.polar(np.radians(sagittarius_angles), sagittarius_radii)
# plt.polar(np.radians(local_angles), local_radii)
# plt.polar(np.radians(perseus_angles), perseus_radii)
# plt.polar(np.radians(outer_angles), outer_radii)


