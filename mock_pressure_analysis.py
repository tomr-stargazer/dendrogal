"""
Automated pressure analyses with fixed distance assumptions.

"""

from __future__ import division

import numpy as np
import astropy
import astropy.units as u

import astrodendro
from demo import orion_demo, perseus_demo, ophiuchus_demo
from assign_physical_values import assign_size_mass_alpha_pressure
from integrated_viewer import IntegratedViewer
from astropy.wcs import wcs
from astrodendro.scatter import Scatter

def fixed_distance_analysis(demo_function, distance, **kwargs):

	d, catalog, datacube_dt_header, metadata = demo_function(**kwargs)

	distance_array = np.ones_like(catalog['_idx']) * u.Quantity(distance, u.pc)

	distance_column = astropy.table.Column(
		data=u.Quantity(distance_array, u.pc),
		name="Distance")
	catalog.add_column(distance_column)

	s, m, v, p = assign_size_mass_alpha_pressure(catalog)

	catalog['size'] = astropy.table.Column(data=s, name='size')
	catalog['mass'] = astropy.table.Column(data=m, name='mass')
	catalog['virial'] = astropy.table.Column(data=v, name='virial')
	catalog['pressure'] = astropy.table.Column(data=p, name='pressure')

	return d, catalog, datacube_dt_header, metadata

def variable_distance_analysis(demo_function, map_idx_distance, **kwargs):

	d, catalog, datacube_dt_header, metadata = demo_function(**kwargs)

	distance_array = np.ones_like(catalog['_idx']) * np.nan

	for key in map_idx_distance.keys():
		idx = key
		distance = map_idx_distance[key]

		trunk = d[idx]
		indices = [structure.idx for structure in trunk.descendants]

		distance_array[np.in1d(catalog['_idx'], indices)] = distance

	distance_column = astropy.table.Column(
		data=u.Quantity(distance_array, u.pc),
		name="Distance")
	catalog.add_column(distance_column)

	s, m, v, p = assign_size_mass_alpha_pressure(catalog)

	catalog['size'] = astropy.table.Column(data=s, name='size')
	catalog['mass'] = astropy.table.Column(data=m, name='mass')
	catalog['virial'] = astropy.table.Column(data=v, name='virial')
	catalog['pressure'] = astropy.table.Column(data=p, name='pressure')

	return d, catalog, datacube_dt_header, metadata

def fixed_perseus_analysis(downsample_factor=1, min_npix=20, **kwargs):
	return fixed_distance_analysis(perseus_demo, 250, downsample_factor=downsample_factor, min_npix=min_npix, **kwargs)

perseus_map = {}
for idx in [64, 296, 373, 499, 392, 63, 11]:
	perseus_map[idx] = 250

def perseus_analysis():
	return variable_distance_analysis(perseus_demo, perseus_map, downsample_factor=1, min_npix=20)

def L1448_analysis():
	return variable_distance_analysis(perseus_demo, {417: 250}, downsample_factor=1, min_npix=20)

alyssa_perseus_map = {966: 250}

def alyssa_perseus_analysis():
	return variable_distance_analysis(perseus_demo, alyssa_perseus_map, downsample_factor=1, min_npix=10)

def ophiuchus_analysis(downsample_factor=1, min_npix=20):
	return fixed_distance_analysis(ophiuchus_demo, 125, downsample_factor=downsample_factor, min_npix=min_npix)

orion_map_idx_distance = {47: 450, 58: 450, 158: 800} # Orion A, Orion B, Monoceros
def orion_analysis():
	return variable_distance_analysis(orion_demo, orion_map_idx_distance, downsample_factor=1, min_npix=20)
