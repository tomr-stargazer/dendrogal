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

    # let's permute. and downsample.
    new_header['naxis'+str(pt[0]+1)] = header['naxis1']
    new_header['naxis'+str(pt[1]+1)] = header['naxis2']
    new_header['naxis'+str(pt[2]+1)] = header['naxis3']

    new_header['ctype'+str(pt[0]+1)] = header['ctype1']
    new_header['ctype'+str(pt[1]+1)] = header['ctype2']
    new_header['ctype'+str(pt[2]+1)] = header['ctype3']

    new_header['crval'+str(pt[0]+1)] = header['crval1']
    new_header['crval'+str(pt[1]+1)] = header['crval2']
    new_header['crval'+str(pt[2]+1)] = header['crval3']

    new_header['cdelt'+str(pt[0]+1)] = header['cdelt1']
    new_header['cdelt'+str(pt[1]+1)] = header['cdelt2']
    new_header['cdelt'+str(pt[2]+1)] = header['cdelt3']

    return new_datacube, new_header


def interpolate_data():
    pass

def moment_mask_data():
    pass

