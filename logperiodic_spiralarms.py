""" 
logperiodic_spiralarms.py

"""

from __future__ import division

import numpy as np

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

    lnR = np.log(Rref) - (beta - betaref)*np.tan(psi)
    spiral_radius = np.exp(lnR)

    return spiral_radius

