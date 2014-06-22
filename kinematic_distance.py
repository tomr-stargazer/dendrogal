"""
This is an attempt at re-implementing a Brand & Blitz kinematic distance tool
in Python, based on Chris Beaumont's kdist.pro and using the new fancy parameters 
from Mark Reid et al. (2014).

"""

from __future__ import division

import numpy as np

import astropy.table
from astropy.coordinates import Galactic
import astropy.units as u
import astropy.constants as c

brand_blitz_parameters = dict()

def make_blitz_distance_column(catalog, nearfar='near'):
    """ Makes a Blitz distance column. """

    nearfar_dict = {'near': 0, 'far': 1}

    distance_list = []

    for row in catalog:
    	blitz_output = brand_blitz_kinematic_distance(row['x_cen'], row['v_cen'])

    	# Sometimes the "near" distance is geometrically negative - obviously, use the far one.
    	if blitz_output[0] < 0:
    		if blitz_output[1] < 0:
    			blitz_output = [np.nan, np.nan]
    		else:
    			blitz_output[0] = blitz_output[1]

    	distance_list.append(blitz_output[nearfar_dict[nearfar]])

    distance_column = astropy.table.Column(
        data=distance_list, 
        name="Distance")

    distance_column.units = u.kpc    
    return distance_column

def assign_blitz_distances(catalog, nearfar='near'):
    # let's pretend these are in... parsecs.
    distance_column = astropy.table.Column(
        data=make_blitz_distance_column(catalog, nearfar), 
        name="Distance")

    distance_column.units = u.kpc
    catalog.add_column(distance_column)

    # misleading to return it, since it's modified in-place
    return catalog




def brand_blitz_kinematic_distance(longitude, velocity, galaxy_parameters=brand_blitz_parameters):
	"""
	Computes a kdist for a given (l,b) pair.

	- notation:
	 v : rotation velocity
	 vo: Solar rotation velocity
	 vr: Radial velocity wrt LSR
	 r : Galactocentric radius
	 ro: Solar galactocentric radius

	rotation curve parameters from Brand and Blitz 1993
	 v/vo = a1 * (R/Ro)^a2 + a3
	 vr = sin(l) * vo * [ a1 * (r/ro)^(a2-1) + a3*(r/ro)^-1 - 1]

	"""
	l = np.radians(longitude) % 2*np.pi

	a1=1
	a2=0
	a3=0

	#- the following are from Reid et al. (2014)
	ro=8.34
	vo=240 # aka theta_0

	r = (np.arange(2000) + 1) / 2000 * 2 * ro
	root = np.sin(l) * vo * (a1 * (r/ro)**(a2-1) + a3 * (r/ro)**(-1) - 1) - velocity

	#find the zero crossing
	root *= np.roll(root,1)
	root[0] = 1
	root[-1] = 1

	hit = np.where(root < 0)

	if len(hit[0]) < 1:
		print 'Cannot determine galactocentric distance'
		print 'returned 0 (NaN) solutions'
		return [np.nan, np.nan]
	r=r[hit[0]][0]

	warnmsg = 'Warning: The galactocentric distance extrapolates the measured galactic roation curve'
	if (r < 2 or r > 17):
		print warnmsg

	# Having determined "r", we now need to solve the triangle and find
	# "d". If "l" is acute, there are (usually) two possible solutions; if
	# "l" is obtuse, then there should be one unique solution.

	# Case I: obtuse l - 1 solution

	if (l >= np.pi/2) and (l < 3*np.pi/2): 

	   # "Gamma" and "Beta" refer to the angles opposite "ro" and "d", 
	   # respectively. (I made a triangle to sketch this all out.)
	   gamma = np.arcsin( ro / r * np.sin(l))
	   beta = np.pi - ( l + gamma)
	   d = r * np.sin(beta) / np.sin( l )

	   print 'returned one solution'
	   return [d, d]

	# Case II: acute l - 2 or 0 solutions

	else:

	   rmin = ro*np.cos(l)
	   dr = np.sqrt(r**2-(ro*np.sin(l))**2)
	   if ~np.isfinite(dr):
	      print 'motion cannot be reproduced via galactic rotation'
	      print 'returned 0 (NaN) solutions'
	      return [np.nan, np.nan]
	   else:
	      print 'returned 2 solutions'
	      return [rmin-dr,rmin+dr]


