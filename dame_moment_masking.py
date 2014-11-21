""" 
This is an implementation of Tom Dame's datacube moment-masking technique. 

A shortcoming of clipping is that it determines whether or not a particular 
spectral channel contains real emission based only on the intensity level in
that channel, whereas a seasoned radio astronomer would, in addition, look 
for spatial and velocity coherence of the signal. Essentially, moment masking
is a refinement of the clipping method in which the coherence of the 
signal---determined from a smoothed version of the data cube---is considered 
in determining which peaks are real (i.e., not blanked).

Described in detail here:
http://www.cfa.harvard.edu/rtdc/CO/MomentMasking/
and here:
http://arxiv.org/abs/1101.1499

"""

from __future__ import division

import numpy as np

from FITS_tools.cube_regrid import gsmooth_cube

def integer_to_tuple(x):
    if x >= 0:
        return (x, 0)
    else:
        return (0, -x)

def integer_to_slice(x):
    if x > 0:
        return slice(None, -x)
    elif x == 0:
        return slice(None, None)
    else:
        return slice(-x, None)

def roll_cube(cube, roll_tuple):

    if len(roll_tuple) != len(cube.shape):
        raise ValueError("`roll_tuple` must match the data cube shape")

    pad_width_tuple = tuple([integer_to_tuple(x) for x in roll_tuple])
    slices = tuple([integer_to_slice(x) for x in roll_tuple])

    shifted_cube = np.pad(cube, pad_width_tuple, mode='constant', constant_values=(np.nan,np.nan))[slices]

    return shifted_cube


def moment_mask(cube, rms_noise, velocity_smoothing=2, spatial_smoothing=2, clip_at_sigma=5):

    # cube : T (v, x, y)

    # T_s (v, x, y)
    smooth_cube = gsmooth_cube(cube, [spatial_smoothing,spatial_smoothing, velocity_smoothing], kernelsize_mult=1)

    smoothed_rms_noise = 1/np.sqrt(spatial_smoothing*spatial_smoothing*velocity_smoothing) * rms_noise

    # T_c
    clipping_level = clip_at_sigma*smoothed_rms_noise

    # M (v, x, y)
    mask_cube = np.zeros_like(cube)

    mask_cube[smooth_cube > clipping_level] = 1

    ns = spatial_smoothing
    nv = velocity_smoothing

    for dx in range(-ns+1, ns):

        for dy in range(-ns+1, ns):

            for dv in range(-nv+1, nv):

                rolled_smooth_cube = roll_cube(smooth_cube, (dv, dx, dy)) # assumes data has been transposed (2, 0, 1)

                mask_cube[rolled_smooth_cube > clipping_level] = 1

    # T_M (v, x, y)
    moment_masked_cube = mask_cube * cube

    return moment_masked_cube
    # , mask_cube, smooth_cube, clipping_level
