"""
This is a demo to open, downsample, and dendrogram some sample data.

The data are allowed to be the "actual" thing we're trying to dendrogram I guess.

"""

from __future__ import division

import numpy as np

import astropy
import astrodendro
import astropy.units as u

from astropy import wcs

from astropy.io.fits import getdata

data_path = "/Users/tsrice/Dropbox/college/Astro99/DATA/"

def cogal_downsampled_demo(downsample_factor=4):

    if type(downsample_factor) != int:
        raise TypeError("Downsample factor must be an Integer")
    elif downsample_factor < 1:
        raise ValueError("Downsample factor must be positive")
    elif downsample_factor > 100:
        raise ValueError("Downsample factor must be less than 100")
    
    df = downsample_factor

    print "loading data: ..."
    cogal, cogal_header = getdata(data_path+'COGAL_all_mom.fits', memmap=True,
                                  header=True)

    cogal_downsampled = cogal[::df, ::df, ::df]

    v_scale = cogal_header['cdelt1'] * df
    v_unit = u.km / u.s
    l_scale = cogal_header['cdelt2'] * df
    b_scale = cogal_header['cdelt3'] * df

    print "computing dendrogram: ..."
    d = astrodendro.Dendrogram.compute(
        cogal_downsampled.transpose(2,0,1), 
        min_value=0.01, min_delta=0.005, min_npix=2000//df**3, verbose=True)

    d.viewer(galactic=True)

    metadata = {}
    metadata['data_unit'] = u.Jy # This is a LIE but Kelvin don't work yet
    metadata['spatial_scale'] = b_scale * u.deg
    metadata['velocity_scale'] = v_scale * v_unit
    metadata['vaxis'] = 1 # because we transposed from (v, l, b) to (2, 0, 1)
                          # for plotting reasons

    catalog = astrodendro.ppv_catalog(d, metadata)

    return d, catalog, cogal_header, metadata
