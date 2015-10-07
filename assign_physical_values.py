"""
This is a "first crack" at assigning physical values to dendrogram parameters.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy
import astrodendro

import astropy.units as u
import astropy.constants as c

try:
    from reid_distance_assigner import make_reid_distance_column
except:
    pass


def assign_distances(catalog, nearfar='near'):
    """ 
    Assigns distances to all entries in the catalog. 

    Currently very wrong -- it's the "stupidest thing you could do".

    """

    # Stupid way to assign distances: random manipulations of x, y, v!

    # let's pretend these are in... parsecs.
    distance_column = astropy.table.Column(
        data=make_reid_distance_column(catalog, nearfar)['D_k'], 
        name="Distance")

    catalog.add_column(distance_column)

    # misleading to return it, since it's modified in-place
    return catalog

def assign_size_mass_alpha_pressure(catalog):
    """
    Assigns physical parameters size, mass, virial parameter, and pressure.

    Does this to a catalog that has distances in it.

    """

    if 'Distance' not in catalog.colnames:
        try:
            catalog['Distance'] = catalog['D_k']
        except Exception:
            raise ValueError("Catalog must have distance!")

    sky_radius = catalog['radius'].data * catalog['radius'].unit
    distance = catalog['Distance'].data * catalog['Distance'].unit
    size = sky_radius.to(u.rad).value * distance
    flux_kelvin_kms_deg2 = catalog['flux_kelvin_kms_deg2'].data * catalog['flux_kelvin_kms_deg2'].unit
    sigma_v = catalog['v_rms'].data * catalog['v_rms'].unit

    luminosity = flux_kelvin_kms_deg2.to(u.K * u.km/u.s * u.steradian) * distance**2 / (1 * u.steradian)

    # X_2 is the ratio between the assumed X_CO value, and the nominal value 2e20 cm-2 (K km s-1)-1 
    X2 = 1.0

    mass = 4.4 * X2 * luminosity.to(u.K * u.km/u.s * u.pc**2).value * u.M_sun

    eta = 1.9

    virial_parameter = 5 * eta * sigma_v**2 * size / (mass * c.G)

    # Pressure is calculated as a kinetic energy density: mass density times velocity squared
    pressure = mass * sigma_v**2 / size**3
    pressure_by_boltzmann = pressure / c.k_B

    return size.to('pc'), mass, virial_parameter.decompose(), pressure_by_boltzmann.to('K cm-3')

def assign_galactocentric_coordinates(catalog, galactic_center_distance=8.340*u.kpc):
    """
    Assigns galactocentric coordinates. 

    Assumes an R_0 = 8340 pc, re: Reid et al. 2013

    This function is derived from the IDL code http://www.astro.virginia.edu/~dln5q/research/idl/lbd2xyz.pro
    by D. Nidever,
    translated into python/astropy

    """

    if 'Distance' not in catalog.colnames:
        try:
            catalog['Distance'] = catalog['D_k']
        except Exception:
            raise ValueError("Catalog must have distance!")

    R_0 = u.Quantity(galactic_center_distance, u.kpc)

    distance = catalog['Distance'].data * catalog['Distance'].unit
    lrad = (catalog['x_cen'] * u.deg).to(u.rad).value
    brad = (catalog['y_cen'] * u.deg).to(u.rad).value

    x = distance * np.sin(0.5*np.pi - brad) *np.cos(lrad) - R_0
    y = distance * np.sin(0.5*np.pi - brad) *np.sin(lrad)
    z = distance * np.cos(0.5*np.pi - brad) 

    return x.to('kpc'), y.to('kpc'), z.to('kpc')

    """
    astropy.table.Column(
        data=catalog['radius'] * (catalog['Distance']*1000) * 3600, 
        name="Size")

    mass = astropy.table.Column(
        data=4.4 * catalog['flux'] / (210099.2) * (catalog['Distance']*1000)**2,
        name="Mass")

    eta = 1.91
    G = 1/232.

    factors = 5 * eta * (catalog['v_rms'])**2 / G
    
    virial_parameter = astropy.table.Column(
        data=factors * size.data / mass.data, 
        name="Virial_Parameter")
    
    catalog.add_columns([size, mass, virial_parameter])

    return catalog
    """

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
