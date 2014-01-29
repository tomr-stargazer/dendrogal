"""
This is a "first crack" at assigning physical values to dendrogram parameters.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy
import astrodendro


def assign_distances(catalog):
    """ 
    Assigns distances to all entries in the catalog. 

    Currently very wrong -- it's the "stupidest thing you could do".

    """

    # Stupid way to assign distances: random manipulations of x, y, v!

    # let's pretend these are in... parsecs.
    distance_column = astropy.table.Column(
        "Distance", 
        60*catalog['v_cen'] + 4*catalog['x_cen'] - 3*catalog['y_cen'])

    catalog.add_column(distance_column)

    # misleading to return it, since it's modified in-place
    return catalog

def assign_size_mass_alpha(catalog):
    """
    Assigns physical parameters size, mass, and virial parameter.

    Does this to a catalog that has distances in it.

    """

    if 'Distance' not in catalog.colnames:
        raise ValueError("Catalog must have distance!")

    size = astropy.table.Column(
        "Size", catalog['radius'] * (catalog['Distance']*1000) * 0.00218)

    mass = astropy.table.Column(
        "Mass", 
        4.4 * catalog['flux'] / (210099.2) * (catalog['Distance']*1000)**2)

    eta = 1.91
    G = 1/232.

    factors = 5 * eta * (catalog['v_rms'])**2 / G
    
    virial_parameter = astropy.table.Column(
        "Virial_Parameter", factors * size.data / mass.data)
    
    catalog.add_columns([size, mass, virial_parameter])

    return catalog

"""
# Don't forget - kdist returns values in kpc, but we want
# parsecs, so we multiply by 1e3.
# AND DON'T FORGET to include our angular resolution (1/8 degree)
# Where in virial does that come in? the "R" value I suspect.

# Don't use these virial estimates until I've worked them out by hand.
# near_virial *= 1./(((catalog.kdist_near)*1d3)^2 * 0.00218) / (210099.2d)
# far_virial *= 1./(((catalog.kdist_far)*1d3)^2 * 0.00218) / (210099.2d)

  eta = 1.91              # correct for concentration of R. see rosolowsky 2008
  G = 1/232.              # in units of km/s and parcsecs. See Solomon 1987

  factors = 5 * eta * (catalog.sig_v)^2 / G

  near_virial = factors * size_near / mass_near
    

# And then masses, which are based on flux and distance
# M = 4.4 * X_2 * L_co
# L_co = F d^2
# where F is in K km s^-1 sr
# and it's tricky to convert from px to steradians. 
# We'd multiply if it were "per" steradian, but steradian is on
# top, so we divide by the number of px elements in a steradian.
# (I hope.)

# flux = catalog.flux / (210099.2d)

# mass_near = 4.4 * flux * (catalog.kdist_near*1d3)^2
# mass_far = 4.4 * flux * (catalog.kdist_far*1d3)^2


# size in pc = size in px * arcseconds per px * AU per arcsecond * pc per AU

# catalog.size_near, eventually
size_near = catalog.sig_r * (catalog.kdist_near*1d3) * 0.00218
catalog.extra1 = size_near
# catalog.size_far, eventually
size_far = catalog.sig_r * (catalog.kdist_far*1d3) * 0.00218
catalog.extra2 = size_far


"""
