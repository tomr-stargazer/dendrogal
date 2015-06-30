"""
Prunes the catalog columns down to what Alyssa wants.

"""

from __future__ import division

import astropy.table
from astropy.coordinates import Galactic
import astropy.units as u


def prune_columns_and_add_radec(catalog):

    new_catalog = astropy.table.Table()

    new_catalog['l'] = catalog['x_cen']
    new_catalog['b'] = catalog['y_cen']
    new_catalog['v'] = catalog['v_cen']

    galactic_coord = Galactic(
        l=catalog['x_cen'], b=catalog['y_cen'], unit=(u.deg, u.deg))

    new_catalog['RA'] = galactic_coord.fk5.ra.deg
    new_catalog['Dec'] = galactic_coord.fk5.dec.deg

    new_catalog['v_rms'] = catalog['v_rms']

    new_catalog['distance'] = catalog['distance']

    new_catalog['error_distance_plus'] = catalog['error_distance_plus']
    new_catalog['error_distance_minus'] = catalog['error_distance_minus']

    new_catalog['mass'] = catalog['mass']

    new_catalog['error_mass_plus'] = catalog['error_mass_plus']
    new_catalog['error_mass_minus'] = catalog['error_mass_minus']

    new_catalog['size'] = catalog['size']

    new_catalog['error_size_plus'] = catalog['error_size_plus']
    new_catalog['error_size_minus'] = catalog['error_size_minus']

    new_catalog['x_gal'] = catalog['x_gal']
    new_catalog['y_gal'] = catalog['y_gal']
    new_catalog['z_gal'] = catalog['z_gal']

    return new_catalog