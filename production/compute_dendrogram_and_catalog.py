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

    Returns
    -------

    """

    if min_value is None or min_delta is None or min_npix is None:
        raise ValueError("Please provide min_value, min_delta, and min_npix!")

    datacube_wcs = wcs.wcs.WCS(header)
    datacube_wcs.wcs.bounds_check(pix2world=False, world2pix=False)

    d = astrodendro.Dendrogram.compute(
        datacube, wcs=datacube_wcs, verbose=verbose,
        min_value=min_value, min_delta=min_delta, min_npix=min_npix)

    # We'll only return d; the data and header are un-altered and should be
    # re-used.
    return d

def compute_catalog(d, header):
    """
    Computes a catalog on a dendrogram.

    Parameters
    ----------
    d : astrodendro.dendrogram.Dendrogram
        Dendrogram on which to compute a catalog
    header : astropy.io.fits.header.Header
        Header corresponding exactly to `d.data`. 
        Used for unit and scaling information to calculate flux and 
        world coordinates.

    Returns
    -------
    catalog : astropy.table.table.Table
        Catalog describing the structures in dendrogram `d`
        using units provided in `header`.
    metadata : dict
        Explicitly lists unit properties used in `catalog`
        
    """

    # rough average of North and South telescopes: 
    # see Dame et al. 2001, p. 794. 
    # "The Northern telescope... has a beamwidth of 8'.4 +/- 0'.2.... 
    # The Southern telescope has nearly the same beamwidth to within 
    # the uncertainties: 8'.8 +/- 0'.2"
    beam_size = 8.5 * u.arcmin 
    frequency = 115.27 * u.GHz # http://www.cv.nrao.edu/php/splat/

    # Remember: by this point, the data ought to be in (v, b, l), with FITS header order opposite that
    if 'vel' not in header['ctype3'].lower():
        raise ValueError("CTYPE3 must be velocity - check that the data were permuted correctly")

    v_scale = header['cdelt3']
    v_unit = u.km / u.s
    l_scale = header['cdelt1']
    b_scale = header['cdelt2']

    metadata = {}
    metadata['data_unit'] = u.K
    metadata['spatial_scale'] = b_scale * u.deg
    metadata['velocity_scale'] = v_scale * v_unit
    metadata['wavelength'] = frequency # formerly: (c.c / frequency).to('mm') but now compute_flux can handle frequency in spectral equivalency
    metadata['beam_major'] = beam_size
    metadata['beam_minor'] = beam_size    
    metadata['vaxis'] = 0 # keep it this way if you think the (post-downsample/transposed) input data is (l, b, v)
    metadata['wcs'] = d.wcs

    catalog = astrodendro.ppv_catalog(d, metadata, verbose=True)

    if catalog['flux'].unit.is_equivalent('Jy'):
        # Workaround because flux is computed wrong (see https://github.com/dendrograms/astrodendro/issues/107)

        flux = u.Quantity(catalog['flux'])
        area_exact = u.Quantity(catalog['area_exact']) #.unit*catalog['area_exact'].data

        # average brightness temperature integrated over area_exact
        flux_kelvin = flux.to('K', equivalencies=u.brightness_temperature(area_exact, frequency))
        # flux integrated over area and velocity
        flux_kelvin_kms_sr = flux_kelvin * metadata['velocity_scale'] * area_exact.to(u.steradian)

        catalog['flux_true'] = flux_kelvin_kms_sr

    return catalog, metadata
