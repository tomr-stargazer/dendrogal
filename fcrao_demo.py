""" Let's try to load the FCRAO data. """

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

from demo import resample_3d_variable

data_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/FCRAO_data/")

# test_load 

def test_load(data_file="grs-35-cube.fits"):

    datacube, datacube_header = getdata(data_path+data_file, memmap=True,
                                        header=True)

    return datacube, datacube_header


spatial_ratio = 20.325203100862584
spectral_ratio = 3.05920142439195

