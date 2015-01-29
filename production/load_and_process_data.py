""" 
This file contains the production code to load & transform datacubes.
"""

from __future__ import division

import pickle
import datetime
import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs
from astropy.io.fits import getdata, getheader
import astropy.io.fits as fits

from .config import data_path

def load_data(filename, data_path=data_path):
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

    return datacube, header

def permute_data_to_standard_order(datacube, header, permute_tuple=(2,0,1)):
    """
    Permutes the axes of a datacube. Default case is from (b,l,v) to (v,b,l).

    I chose this permutation order because of how `plt.imshow` handles 2D arrays.

    Parameters
    ----------
    datacube : np.ndarray

    header : astropy.io.fits.header.Header

    permute_tuple : 3-element tuple, default: (2,0,1)
        A tuple expressing how to permute the axes, 
        according to np.ndarray.transpose(). 

    Returns
    -------
    new_datacube : np.ndarray
        The result of running `datacube.transpose(permute_tuple)`.
    new_header : astropy.io.fits.header.Header
        The header with the axes' properties permuted.

    """

    pt = permute_tuple

    new_datacube = datacube.transpose(pt)

    new_header = header.copy()

    # let's permute.
    for keyword in ['naxis', 'ctype', 'crval', 'cdelt', 'crpix']:
        for i in range(3):
            # replace each keyword's value at `i+1` with its value at `pt[i]+1`.
            new_header['{0}{1}'.format(keyword, pt[i]+1)] = header['{0}{1}'.format(keyword, i+1)]

    return new_datacube, new_header


def interpolate_data():
    pass

def moment_mask_data():
    pass

