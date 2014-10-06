""" 
This is an implementation of Tom Dame's datacube interpolation technique. 

"In each spectrum, <= 2 missing channels are filled by linear interpolation. 
In each spatial plane, single missing pixels are filled by linear 
interpolation, first in l direction, then b."

http://www.cfa.harvard.edu/rtdc/CO/CompositeSurveys/

"""

from __future__ import division

import numpy as np

def shift_left(array, n=1):

	return np.pad(array, (0,n), mode='constant', constant_values=(0,np.nan))[n:]

def shift_right(array, n=1):

	return np.pad(array, (n,0), mode='constant', constant_values=(np.nan,0))[:-n]


def interpolate_single(array):
	""" Linearly interpolates single `nan` values based on their neighbors. """

	new_array = np.copy(array)

	# 1. find single nans in the array -- i.e. nans whose neighbors are both NOT nans
	nan_positions = np.where(np.isnan(array) & ~np.isnan(shift_right(array)) & ~np.isnan(shift_left(array)))

	for i in nan_positions[0]:
		# interpolate between neighbors
		new_array[i] = np.mean((new_array[i+1], new_array[i-1]))

	return new_array

def interpolate_double(array):
	""" Linearly interpolates double `nan` values based on their neighbors. Ignores single nans! """

	new_array = np.copy(array)

	# 1. find double nans in the array -- i.e. two nans whose neighbors are both NOT nans.
	# We only count the bottom of the two, to simplify things. (symmetrical so that's good.)
	nan_positions = np.where(np.isnan(array) & np.isnan(shift_left(array,1)) & ~np.isnan(shift_right(array,1)) & ~np.isnan(shift_left(array,2)))

	for i in nan_positions[0]:
		# interpolate between neighbors
		new_array[i] = (2/3 * new_array[i-1]) + (1/3 * new_array[i+2])
		new_array[i+1] = (1/3 * new_array[i-1]) + (2/3 * new_array[i+2])

	return new_array

def interpolate_spectrum(array):
	""" Interpolates a spectrum, first doing single-nans, then double-nans """
	intermediate_array = interpolate_single(array)
	return interpolate_double(intermediate_array)


def interpolate_datacube(data, spectrum_axis=2, lon_axis=1, lat_axis=0):
	""" 
	Interpolate missing values in a datacube according to Dame's prescription.

	Default: assumes (b, l, v) cube.

	"""

	if len(data.shape) != 3:
		raise ValueError("`data` must have dimensions=3")

	if (spectrum_axis != 2) or (lon_axis != 1) or (lat_axis != 0):
		raise NotImplementedError("Function can't yet handle different axis permutations")

	new_data = np.copy(data)

	# interpolate each spectrum, fix in-place
	for l_i in range(new_data.shape[lon_axis]):
		for b_i in range(new_data.shape[lat_axis]):

			old_spectrum = new_data[b_i, l_i, :] # In principle, the order of [b, l, :] should be dependent on spectrum_axis and the like. But I don't know how to do that.
			new_data[b_i, l_i, :] = interpolate_spectrum(old_spectrum)

	# interpolate each longitude slice
	for v_i in range(new_data.shape[spectrum_axis]):
		for b_i in range(new_data.shape[lat_axis]):

			old_lon_slice = new_data[b_i, :, v_i]
			new_data[b_i, :, v_i] = interpolate_single(old_lon_slice)

	# interpolate each latitude slice
	for v_i in range(new_data.shape[spectrum_axis]):
		for l_i in range(new_data.shape[lon_axis]):

			old_lat_slice = new_data[:, l_i, v_i]
			new_data[:, l_i, v_i] = interpolate_single(old_lat_slice)

	return new_data




