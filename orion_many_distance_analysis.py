""" Let's place Orion at lots of distances and see how it looks!

Downsample spatially but not spectrally. Later: introduce noise.

"""

from __future__ import division

import pickle
import os.path
import datetime
import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs
from astropy.io.fits import getdata, getheader
import astropy.io.fits as fits

from FITS_tools.cube_regrid import gsmooth_cube

from demo import resample_3d_variable, downsampled_demo, downsample_and_transpose_data_and_header
from dame_moment_masking import moment_mask

data_path = os.path.expanduser("~/Dropbox/College/Astro99/DATA/")

data_file = 'DHT27_Orion_interp.fits'

def test_orion(data_file=data_file, distance_multiplier=1, **kwargs):

    df = (distance_multiplier, distance_multiplier, 1)

    # should invoke resample_3d_variable
    d, catalog, header, metadata = downsampled_demo(data_file=data_file, data_path=data_path, 
        downsample_factor=df, min_npix=75, min_value=0.05, min_delta=0.1, recenter=False, **kwargs)


    return d, catalog, header, metadata

# then you downsample in various ways, add noise at each "distance", measure what you get out!
# using the dendrogram extraction of course (:

def determine_number_of_structures_as_a_function_of_distance():

    dm_array = np.arange(19)+2

    n_structs = []

    max_v_rms = []

    for i in dm_array:

        d, catalog, header, metadata = test_orion(distance_multiplier=i)

        n_structs.append(len(catalog))

        max_v_rms.append(catalog['v_rms'][catalog['radius'] ==max(catalog['radius'])][0])

    return dm_array, n_structs, max_v_rms



def downsample_addnoise_and_momentmask(datacube, header, downsample_factor, rms_noise):

    smooth_cube = gsmooth_cube(datacube, [downsample_factor, downsample_factor, 1], kernelsize_mult=1)

    downsample_tuple = (downsample_factor, downsample_factor, 1)

    downsampled_datacube, downsampled_header = downsample_and_transpose_data_and_header(smooth_cube, header, downsample_factor=downsample_tuple, resample=True)

    noise_cube = np.random.standard_normal(downsampled_datacube.shape) * rms_noise

    noisy_cube = downsampled_datacube + noise_cube

    moment_masked_cube = moment_mask(noisy_cube, rms_noise, velocity_smoothing=4)

    return moment_masked_cube, downsampled_header







