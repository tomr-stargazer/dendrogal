""" 
This file contains the production code to load & transform datacubes.
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

data_path = os.path.expanduser("~/Dropbox/College/Astro99/DATA/")

def load_data(filename):
    """
    Loads a datacube (and its header) from disk using `astropy.io.fits.getdata`.

    Assumes the data is in `data_path`. Uses memmap=True.

    Parameters
    ----------
    filename : str
        Name of the file within `data_path`.
        ex: "DHT17_Quad2_bw_mom.fits"

    Returns
    -------
    datacube : np.ndarray
        The datacube inside `filename`. 
        DHT datacubes typically load in (b,l,v) order when using numpy.
    header : astropy.io.fits.header.Header

    """

    # memmap is set to True because large files can otherwise slow us down
    datacube, header = getdata(data_path+filename, header=True, memmap=True)

def transpose_data_to_standard_order(datacube, header):
    return datacube, header



def interpolate_data():
    pass

def moment_mask_data():
    pass

