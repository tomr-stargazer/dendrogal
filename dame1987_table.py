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
	'Cepheus': {'l_min':100,
	            'l_max':120,
	            'b_min':11,
	            'b_max':22,
	            'distance':450,
	            'mass':1.9},	              
	'Taurus':  {'l_min':163,
	            'l_max':178,
	            'b_min':-22,
	            'b_max':-9.5,
	            'distance':140,
	            'mass':0.3},
	'Per OB2 1':{'l_min':154,
	             'l_max':162.5,
	             'b_min':-25,
	             'b_max':-7,
	             'distance':350,
	             'mass':1.3},	   
	'Per OB2 2':{'l_min':163,
	             'l_max':171,
	             'b_min':-9,
	             'b_max':-6,
	             'distance':350,
	             'mass':1.3},	   	             
	'Mon OB1': {'l_min':197.5,
	            'l_max':205,
	            'b_min':-1,
	            'b_max':4,
	            'distance':800,
	            'mass':1.6},	   
	'Orion A': {'l_min':208.5,
	            'l_max':218,
	            'b_min':-21,
	            'b_max':-14.5,
	            'distance':500,
	            'mass':1.6},	   
	'Orion B': {'l_min':202.5,
	            'l_max':208,
	            'b_min':-21,
	            'b_max':-6,
	            'distance':500,
	            'mass':1.7},	   
	'Mon R2': {'l_min':210,
	           'l_max':218,
	           'b_min':-14,
	           'b_max':-10,
	           'distance':830,
	           'mass':1.2},   
	'Vela Sheet': {'l_min':272,
	               'l_max':279,
	               'b_min':-3,
	               'b_max':8,
	               'distance':425,
	               'mass':0.8},	   
	'Cham': {'l_min':295,
	         'l_max':305,
	         'b_min':-20,
	         'b_max':-12,
	         'distance':215,
	         'mass':0.1},
	'Coalsack': {'l_min':300,
	             'l_max':307,
	             'b_min':-4,
	             'b_max':3,
	             'distance':175,
	             'mass':0.04},	   
	'G317-4': {'l_min':315,
	           'l_max':320,
	           'b_min':-6,
	           'b_max':-2,
	           'distance':170,
	           'mass':0.03},	   
	'Lupus': {'l_min':333,
	          'l_max':346,
	          'b_min':4,
	          'b_max':22,
	          'distance':170,
	          'mass':0.3},	   
	'rho Oph 1': {'l_min':350,
	              'l_max':2,
	              'b_min':13,
	              'b_max':24,
	              'distance':165,
	              'mass':0.3},	   
	'rho Oph 2': {'l_min':356,
	              'l_max':5,
	              'b_min':3,
	              'b_max':12.5,
	              'distance':165,
	              'mass':0.3},	   
	'R CrA': {'l_min':357,
	          'l_max':4,
	          'b_min':-22,
	          'b_max':-14,
	          'distance':150,
	          'mass':0.03},	   
}

for key in table2:
	table2[key]['mass'] *= 1e5

def convert_dame_table_to_standard_form(table):

	for key in table:

		row = table[key]

		if (row['l_min'] > 270) and (row['l_max'] < 90):
			row['l_max'] += 360

		radius = np.sqrt((row['l_max']-row['l_min'])*
			             (row['b_max']-row['b_min']) )

		table[key]['radius'] = radius * u.deg

		l_cen = (row['l_max'] + row['l_min'])/2 % 360
		b_cen = (row['b_max'] + row['b_min'])/2		

		table[key]['l'] = l_cen * u.deg
		table[key]['b'] = b_cen * u.deg

	return table