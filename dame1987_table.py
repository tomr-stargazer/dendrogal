"""
This is a hand-copied version of the Dame 1987 table of molecular clouds

"""

from __future__ import division

import numpy as np

import astropy.units as u

table2 = {
	'Aquila Rift 1': {'l_min':18.5,
	                  'l_max':34,
	                  'b_min':-6,
	                  'b_max':10,
	                  'distance':200,
	                  'mass':1.5},
	'Aquila Rift 2': {'l_min':34,
	                  'l_max':44,
	                  'b_min':-4,
	                  'b_max':4,
	                  'distance':200,
	                  'mass':1.5},
	'Cloud A': {'l_min':44,
	            'l_max':49.5,
	            'b_min':-4,
	            'b_max':2,
	            'distance':500,
	            'mass':0.4},
	'Cloud B': {'l_min':44,
	            'l_max':54,
	            'b_min':-4,
	            'b_max':5,
	            'distance':300,
	            'mass':0.4},
	'Cloud C': {'l_min':50,
	            'l_max':55,
	            'b_min':-1,
	            'b_max':3.5,
	            'distance':500,
	            'mass':0.3},
	'Vul Rift': {'l_min':54,
	             'l_max':63,
	             'b_min':-3,
	             'b_max':5,
	             'distance':400,
	             'mass':0.8},
	'Cyg Rift': {'l_min':63,
	             'l_max':86.5,
	             'b_min':-4,
	             'b_max':4,
	             'distance':700,
	             'mass':8.6},
	'Cyg OB7':  {'l_min':87,
	             'l_max':99,
	             'b_min':-3,
	             'b_max':8,
	             'distance':800,
	             'mass':7.5},
	'Lindblad Ring': {'l_min':100,
	                  'l_max':164,
	                  'b_min':-4,
	                  'b_max':10,
	                  'distance':300,
	                  'mass':1.6},
	'"-12 km/s"':  {'l_min':102,
	                'l_max':161,
	                'b_min':-4,
	                'b_max':10,
	                'distance':800,
	                'mass':8.7},
	'Taurus':  {'l_min':163,
	            'l_max':178,
	            'b_min':-22,
	            'b_max':-9.5,
	            'distance':140,
	            'mass':0.3},
}

def convert_dame_table_to_standard_form(table):

	for key in table:

		row = table[key]

		radius = np.sqrt((row['l_max']-row['l_min'])*
			             (row['b_max']-row['b_min']) )

		table[key]['radius'] = radius * u.deg

		l_cen = (row['l_max'] + row['l_min'])/2
		b_cen = (row['b_max'] + row['b_min'])/2		

		table[key]['l'] = l_cen * u.deg
		table[key]['b'] = b_cen * u.deg

	return table