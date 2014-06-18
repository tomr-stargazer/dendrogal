"""
A module to help assign distances to objects of known coordinates and size

We're setting up a correspondence between objects of known coordinates, size and distance,
and dendrogram objects.

"""

from __future__ import division

import numpy as np

from astropy import units as u
from astropy import constants as c

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

def assign_distance_thingy(lookup, catalog):
	""" Takes a dict of stuff and assigns distances. """

	# for each item in thing: pick a subset of the catalog closest
	for key in lookup:

		nearby_indices = ((np.abs(catalog['x_cen'] - lookup[key]['l'].value) < lookup[key]['radius'].value) &
		                  (np.abs(catalog['y_cen'] - lookup[key]['b'].value) < lookup[key]['radius'].value) )

		size_difference = np.abs(catalog['radius'] - lookup[key]['radius'])

		match = catalog[nearby_indices][ size_difference[nearby_indices] == np.min(size_difference[nearby_indices]) ]

		print "Match! name: {0} to idx: {1}".format(key, match['_idx'].data[0])

		# Eventually our goal is to assign that DISTANCE to the match, but let's not get ahead of ourselves



