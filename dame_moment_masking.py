""" 
This is an implementation of Tom Dame's datacube moment-masking technique. 

As Tom Dame says:
"A shortcoming of clipping is that it determines whether or not a particular 
spectral channel contains real emission based only on the intensity level in
that channel, whereas a seasoned radio astronomer would, in addition, look 
for spatial and velocity coherence of the signal. Essentially, moment masking
is a refinement of the clipping method in which the coherence of the 
signal---determined from a smoothed version of the data cube---is considered 
in determining which peaks are real (i.e., not blanked)."

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


def moment_mask(cube, rms_noise, smoothed_rms_noise=None, velocity_smoothing=2, spatial_smoothing=2, clip_at_sigma=5, velocity_axis=0):
    """
    Moment-masks a cube according to Dame 2011 prescription.

    Parameters
    ----------
    cube : numpy.ndarray
        PPV datacube with velocity axis first (i.e. numpy zero-th axis).
        Dame (2011) calls this "T (v,x,y)".
    rms_noise : float
        RMS noise per channel in `cube`. 
        Only used if `smoothed_rms_noise` is not given.
    smoothed_rms_noise : float or None, optional
        I recommend you provide this value.
        If provided, then this is the rms noise in the smooth cube which
        is used (along with `clip_at_sigma`) for clipping the mask cube.
    velocity_smoothing : float, optional
        The data will be convolved along the velocity axis with a Gaussian
        that has a FWHM of `velocity_smoothing` times the data's velocity 
        resolution.
    spatial_smoothing : float, optional
        The data will be convolved along each spatial axis with a Gaussian
        that has a FWHM of `spatial_smoothing` times the data's angular 
        resolution.
    clip_at_sigma : float, optional
        What factor times `smoothed_rms_noise` should the mask cube clip at?
        Dame (2011) recommends that `clip_at_sigma`=5.

    Returns
    -------
    moment_masked_cube : numpy.ndarray
        A noise-suppressed version of `cube`.
        Dame (2011) calls this "T_M (v,x,y)".

    """

    if velocity_axis != 0:
        raise NotImplementedError("Velocity axis has to be zero right now.")

    # cube : T (v, x, y)

    # convert between FWHM and Gaussian "sigma"
    velocity_sigma = velocity_smoothing/(2*np.sqrt(2*np.log(2)))
    spatial_sigma = spatial_smoothing/(2*np.sqrt(2*np.log(2)))

    kernelsize = 3 * [spatial_sigma]
    kernelsize[velocity_axis] = velocity_sigma

    # T_s (v, x, y)
    smooth_cube = gsmooth_cube(cube, kernelsize, kernelsize_mult=4)

    if smoothed_rms_noise is None:
        # yes, square root not cube root: noise goes down as (# pixels binned) squared.
        smoothed_rms_noise = 1/np.sqrt(spatial_smoothing*spatial_smoothing*velocity_smoothing) * rms_noise

    # T_c
    clipping_level = clip_at_sigma*smoothed_rms_noise

    # M (v, x, y)
    mask_cube = np.zeros_like(cube)

    mask_cube[smooth_cube > clipping_level] = 1

    # This nomenclature departs from the Dame (2011) paper: 
    #  My `ns` and `nv` are double that of Dame's.
    #  I do this to simplify the following code.
    ns = spatial_smoothing
    nv = velocity_smoothing

    for dx in range(-ns+1, ns):

        for dy in range(-ns+1, ns):

            for dv in range(-nv+1, nv):

                # The following line needs to be fixed before we can lift the NotImplementedError exception.
                rolled_smooth_cube = roll_cube(smooth_cube, (dv, dx, dy)) # assumes data has been transposed (2, 0, 1)

                mask_cube[rolled_smooth_cube > clipping_level] = 1

    # T_M (v, x, y)
    moment_masked_cube = mask_cube * cube

    return moment_masked_cube
    # , mask_cube, smooth_cube, clipping_level
