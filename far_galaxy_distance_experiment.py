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


def far_galaxy_distance_demo(resample=2, **kwargs):

	d, catalog, header, metadata = demo.cogal_deep_resampled_demo(resample=resample, min_npix=50, **kwargs)

	blitz = make_blitz_distance_column(catalog)
	catalog['Distance'] = blitz

	x, y, z = assign_galactocentric_coordinates(catalog, galactic_center_distance=0)

	catalog['x_galactocentric'] = x
	catalog['y_galactocentric'] = y
	catalog['z_galactocentric'] = z

	s, m, v, p = assign_size_mass_alpha_pressure(catalog)

	catalog['size'] = astropy.table.Column(data=s, name='size')
	catalog['mass'] = astropy.table.Column(data=m, name='mass')
	catalog['virial'] = astropy.table.Column(data=v, name='virial')
	catalog['pressure'] = astropy.table.Column(data=p, name='pressure')	

	print("\n# Run these commands:\n"
		"# (these assume you have called this function like following: )\n"
		"# d, catalog, x, y = near_galaxy_distance_demo(resample=2) \n"
		"dv = d.viewer()\n"
		"iv = IntegratedViewer(d, dv.hub, wcs=y['wcs'].sub([wcs.WCSSUB_CELESTIAL]), cmap='gray_r')\n"
		"dsd = Scatter(d, dv.hub, catalog, 'y_galactocentric', 'x_galactocentric')")

	return d, catalog, header, metadata

