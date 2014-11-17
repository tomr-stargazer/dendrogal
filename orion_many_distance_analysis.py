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

from demo import resample_3d_variable, downsampled_demo

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

    dm_array = np.arange(20)+1

    n_structs = []

    for i in dm_array:

        d, catalog, header, metadata = test_orion(distance_multiplier=i)

        n_structs.append(len(catalog))

    return dm_array, n_structs

