"""
Functions to turn data into dendrograms & catalogs.

"""

from __future__ import division

import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs


def compute_dendrogram(datacube, header, verbose=True,
                       min_value=None, min_delta=None, min_npix=None):
    """
    Computes a dendrogram on input data/header.

    Parameters
    ----------
    datacube : np.ndarray
    header : astropy.io.fits.header.Header
        Used for WCS information (required!)
    min_value : float 
        the minimum value to consider in the dataset - any value lower 
        than this will not be considered in the dendrogram. If you are 
        working with observations, it is likely that you will want to 
        set this to the "detection level," for example 3- or 5-sigma, 
        so that only significant values are included in the dendrogram. 
        By default, all values are used.
    min_delta : float 
        how significant a leaf has to be in order to be considered an 
        independent entity. The significance is measured from the 
        difference between its peak flux and the value at which it is 
        being merged into the tree. If you are working with observational
        data, then you could set this to, e.g., 1-sigma, which means that 
        any leaf that is locally less than 1-sigma tall is combined with 
        its neighboring leaf or branch and is no longer considered a 
        separate entity.
    min_npix : int 
        the minimum number of pixels/values needed for a leaf to be 
        considered an independent entity. When the dendrogram is being 
        computed, and when a leaf is about to be joined onto a branch 
        or another leaf, if the leaf has fewer than this number of pixels, 
        then it is combined with the branch or leaf it is being merged 
        with and is no longer considered a separate entity. By default, 
        this parameter is set to zero, so there is no minimum number of 
        pixels required for leaves to remain independent entities.
    verbose : bool, default True
        Provide a progress bar? If False, the computation will be 
        silent, but for large dendrograms, it can be useful to have an 
        idea of how long the computation will take.

    """

    if min_value is None or min_delta is None or min_npix is None:
        raise ValueError("Please provide min_value, min_delta, and min_npix!")

    datacube_wcs = wcs.wcs.WCS(header)
    datacube_wcs.wcs.bounds_check(pix2world=False)

    d = astrodendro.Dendrogram.compute(
        datacube, wcs=datacube_wcs, verbose=verbose,
        min_value=min_value, min_delta=min_delta, min_npix=min_npix)

    # We'll only return d; the data and header are un-altered and should be
    # re-used.
    return d

def compute_catalog(d, header):

    pass