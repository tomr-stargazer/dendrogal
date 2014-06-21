"""
A module to help assign distances to objects of known coordinates and size

We're setting up a correspondence between objects of known coordinates, size and distance,
and dendrogram objects.

"""

from __future__ import division

import numpy as np

import astropy
from astropy import units as u
from astropy import constants as c

from dame1987_table import convert_dame_table_to_standard_form, table2
import demo
from assign_physical_values import assign_galactocentric_coordinates, assign_size_mass_alpha_pressure
from astrodendro_analysis.integrated_viewer import IntegratedViewer
from astropy.wcs import wcs
from astrodendro.scatter import Scatter

# Base case: you know the coordinate, the size, the distance, and you want to match a dendrogram object to it.

test_lookup = {'Orion': {'l':210.78199346*u.deg, 
	                     'b':-19.26107202*u.deg, 
	                     'radius':1.17093033*u.deg, 
	                     'distance':450*u.pc},
 	           'Perseus': {'l':159.06126878*u.deg,
	                       'b':-20.26022578*u.deg,
	                       'radius':0.99591075*u.deg,
	                       'distance':250*u.pc} 
	           }

def assign_local_distance(lookup, catalog, reset=True):
	""" Takes a dict of local cloud information and assigns distances. """

	if 'Distance' not in catalog.colnames:
		distance_column = astropy.table.Column(data=np.nan*np.ones(len(catalog))*u.pc, name='Distance')  # shape=(2,)
		catalog.add_column(distance_column)
	if reset:
		catalog['Distance'] *= np.nan

	# for each item in thing: pick a subset of the catalog closest
	for key in lookup:

		nearby_indices = ((np.abs(catalog['x_cen'] - lookup[key]['l'].value) < lookup[key]['radius'].value/2) &
		                  (np.abs(catalog['y_cen'] - lookup[key]['b'].value) < lookup[key]['radius'].value/2) )

		size_difference = np.abs(catalog['radius'] - lookup[key]['radius'])

		valid_aspect_ratio = catalog['major_sigma'] / catalog['minor_sigma'] < 10

		valid_indices = nearby_indices & valid_aspect_ratio

		try:
			match_index = (size_difference == np.min(size_difference[valid_indices])) & valid_indices

			print "Match! name: {0} to idx: {1}".format(key, catalog[match_index]['_idx'].data[0])

			catalog['Distance'][match_index] = lookup[key]['distance']

			print "Matched distance: {0}".format(catalog['Distance'][match_index].data)
		except ValueError, e:
			print "No match! name: {0} had {1}".format(key, e)




