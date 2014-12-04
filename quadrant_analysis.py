"""
Parallel to near_galaxy_distance_experiment

Uses kdist

"""

from __future__ import division

import numpy as np

import astropy
from astropy import units as u
from astropy import constants as c

import demo
from reid_distance_assigner import make_reid_distance_column
from kinematic_distance import make_blitz_distance_column
from assign_physical_values import assign_galactocentric_coordinates, assign_size_mass_alpha_pressure
from astrodendro_analysis.integrated_viewer import IntegratedViewer
from astropy.wcs import wcs
from astrodendro.scatter import Scatter


def secondquad_distance_demo(downsample=1, distance='reid', nearfar='near', min_npix=20, min_value=0.05, min_delta=0.3, **kwargs):

    d, catalog, header, metadata = demo.secondquad_demo_mominterp(downsample_factor=downsample, 
        min_npix=min_npix, min_value=min_value, min_delta=min_delta, recenter=False, **kwargs)

    if distance != 'reid':
        blitz = make_blitz_distance_column(catalog)
        catalog['Distance'] = blitz
    else:
        reid_distance = make_reid_distance_column(catalog, nearfar=nearfar)
        catalog['Distance'] = reid_distance['D_k']

    x, y, z = assign_galactocentric_coordinates(catalog, galactic_center_distance=0)

    catalog['x_galactocentric'] = x
    catalog['y_galactocentric'] = y
    catalog['z_galactocentric'] = z

    s, m, v, p = assign_size_mass_alpha_pressure(catalog)

    catalog['size'] = s
    catalog['mass'] = m
    catalog['virial'] = v
    catalog['pressure'] = p

    # disqualify!

    disqualified = ((catalog['mass'] > 1e8) | 
                    (catalog['mass'] < 5e3) | 
                    (catalog['major_sigma'] > 10) | 
                    (catalog['v_rms'] > 30) | 
                    (catalog['size'] > 1000) |
                    (np.abs(catalog['v_cen']) < 13) |
                    (np.abs(catalog['x_cen'] - 180) < 10) |
                    (catalog['area_exact'] > 50) )

    catalog['Distance'][disqualified] = np.nan
    catalog['v_rms'][disqualified] = np.nan
    catalog['size'][disqualified] = np.nan
    catalog['virial'][disqualified] = np.nan
    catalog['mass'][disqualified] = np.nan


    print("\n# Run these commands:\n"
        "# (these assume you have called this function like following: )\n"
        "# d, catalog, x, y = near_galaxy_distance_demo(resample=2) \n"
        "dv = d.viewer()\n"
        "iv = IntegratedViewer(d, dv.hub, wcs=y['wcs'].sub([wcs.WCSSUB_CELESTIAL]), cmap='gray_r')\n"
        "dsd = Scatter(d, dv.hub, catalog, 'y_galactocentric', 'x_galactocentric')")

    return d, catalog, header, metadata


def thirdquad_distance_demo(downsample=1, distance='reid', nearfar='near', min_npix=20, min_value=0.05, min_delta=0.3, **kwargs):

    d, catalog, header, metadata = demo.thirdquad_demo_mominterp(downsample_factor=downsample, 
        min_npix=min_npix, min_value=min_value, min_delta=min_delta, recenter=False, **kwargs)

    if distance != 'reid':
        blitz = make_blitz_distance_column(catalog)
        catalog['Distance'] = blitz
    else:
        reid_distance = make_reid_distance_column(catalog, nearfar=nearfar)
        catalog['Distance'] = reid_distance['D_k']

    x, y, z = assign_galactocentric_coordinates(catalog, galactic_center_distance=0)

    catalog['x_galactocentric'] = x
    catalog['y_galactocentric'] = y
    catalog['z_galactocentric'] = z

    s, m, v, p = assign_size_mass_alpha_pressure(catalog)

    catalog['size'] = s
    catalog['mass'] = m
    catalog['virial'] = v
    catalog['pressure'] = p

    # disqualify!

    disqualified = ((catalog['mass'] > 1e8) | 
                    (catalog['mass'] < 5e3) | 
                    # (catalog['major_sigma'] > 10) | 
                    # (catalog['v_rms'] > 30) | 
                    # (catalog['size'] > 1000) |
                    # (np.abs(catalog['v_cen']) < 13) |
                    # (np.abs(catalog['x_cen'] - 180) < 10) |
                    (catalog['area_exact'] > 50) )

    catalog['Distance'][disqualified] = np.nan
    catalog['v_rms'][disqualified] = np.nan
    catalog['size'][disqualified] = np.nan
    catalog['virial'][disqualified] = np.nan
    catalog['mass'][disqualified] = np.nan


    print("\n# Run these commands:\n"
        "# (these assume you have called this function like following: )\n"
        "# d, catalog, x, y = near_galaxy_distance_demo(resample=2) \n"
        "dv = d.viewer()\n"
        "iv = IntegratedViewer(d, dv.hub, wcs=y['wcs'].sub([wcs.WCSSUB_CELESTIAL]), cmap='gray_r')\n"
        "dsd = Scatter(d, dv.hub, catalog, 'y_galactocentric', 'x_galactocentric')")

    return d, catalog, header, metadata


def firstquad_distance_demo(downsample=1, distance='reid', min_npix=20, min_value=0.05, min_delta=0.3, **kwargs):

    d, catalog, header, metadata = demo.firstquad_demo_mominterp(downsample_factor=downsample, 
        min_npix=min_npix, min_value=min_value, min_delta=min_delta, recenter=False, **kwargs)

    near_distance_table = make_reid_distance_column(catalog, nearfar='near')
    far_distance_table = make_reid_distance_column(catalog, nearfar='far')

    near_distance_column = near_distance_table['D_k']
    far_distance_column = far_distance_table['D_k']

    best_distance = np.zeros_like(near_distance_table['D_k'])

    sky_radius = u.Quantity(catalog['radius'].data * catalog['radius'].unit)
    near_distance = u.Quantity(near_distance_column)
    near_size = sky_radius.to(u.rad).value * near_distance

    far_distance = u.Quantity(far_distance_column)
    far_size = sky_radius.to(u.rad).value * far_distance

    quad2_fit_constant = 0.48293812090592952
    quad2_fit_power = 0.56796770148326814

    expected_size =  (1/quad2_fit_constant * catalog['v_rms'].data)**(1/quad2_fit_power) * u.pc

    use_near_distance = (np.abs(near_size - expected_size) <= np.abs(far_size - expected_size))
    use_far_distance = (np.abs(near_size - expected_size) > np.abs(far_size - expected_size))

    best_distance[use_near_distance] = near_distance[use_near_distance]
    best_distance[use_far_distance] = far_distance[use_far_distance]



    catalog['Distance'] = best_distance

    x, y, z = assign_galactocentric_coordinates(catalog, galactic_center_distance=0)

    catalog['x_galactocentric'] = x
    catalog['y_galactocentric'] = y
    catalog['z_galactocentric'] = z

    s, m, v, p = assign_size_mass_alpha_pressure(catalog)

    catalog['size'] = s
    catalog['mass'] = m
    catalog['virial'] = v
    catalog['pressure'] = p

    # disqualify!

    disqualified = ((catalog['mass'] > 1e7) | 
                    (catalog['mass'] < 5e3) | 
                    (catalog['major_sigma'] > 10) | 
                    (catalog['v_rms'] > 30) | 
                    (catalog['size'] > 1000) |
                    # (np.abs(catalog['v_cen']) < 13) |
                    # (np.abs(catalog['x_cen'] - 180) < 10) |
                    (catalog['area_exact'] > 50) )

    catalog['Distance'][disqualified] = np.nan
    catalog['v_rms'][disqualified] = np.nan
    catalog['size'][disqualified] = np.nan
    catalog['virial'][disqualified] = np.nan
    catalog['mass'][disqualified] = np.nan


    print("\n# Run these commands:\n"
        "# (these assume you have called this function like following: )\n"
        "# d, catalog, x, y = near_galaxy_distance_demo(resample=2) \n"
        "dv = d.viewer()\n"
        "iv = IntegratedViewer(d, dv.hub, wcs=y['wcs'].sub([wcs.WCSSUB_CELESTIAL]), cmap='gray_r')\n"
        "dsd = Scatter(d, dv.hub, catalog, 'y_galactocentric', 'x_galactocentric')")

    return d, catalog, header, metadata, near_distance_table, far_distance_table

