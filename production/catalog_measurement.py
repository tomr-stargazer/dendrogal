"""
Functions to measure overall cloud statistics, given a sanitized input cloud catalog.

"""

from __future__ import division

import numpy as np

from scipy.odr import RealData, Model, ODR

def powerlaw_function(B, x):
    return B[0] * (x**B[1])

def size_linewidth_slope(catalog):
	"""
	Measures a catalog's slope of sizes and linewidths using scipy.odr.

	ODR means "orthogonal distance regression", and is a way of fitting
	models to data where both the x and y values have scatter.

	Parameters
	----------
	catalog : astropy.table.table.Table
	    Table containing columns for 'size' and 'v_rms'.

	Returns
	-------
	odr_output : scipy.odr.odrpack.Output
	    Output of the scipy ODR routine. The fit parameters
	    are stored in odr_output.beta .

	"""

	if 'size' not in catalog.colnames or 'v_rms' not in catalog.colnames:
		raise ValueError("'size' and 'v_rms' must be columns in `catalog`!")

	size_linewidth_data = RealData(catalog['size'], catalog['v_rms'])

	powerlaw_model = Model(powerlaw_function)
	# The provided `beta0` is a throwaway guess that shouldn't change the outcome.
	odr_object = ODR(size_linewidth_data, powerlaw_model, beta0=[1.,1.])
	odr_output = odr_object.run()

	return odr_output

def cumulative_massfunction_slope():
	pass






