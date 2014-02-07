"""
This is a demo to open, downsample, and dendrogram some sample data.

The data are allowed to be the "actual" thing we're trying to dendrogram I guess.

"""

from __future__ import division

import numpy as np

import astropy
import astrodendro
import astropy.units as u

from astropy.io.fits import getdata

data_path = "/Users/tsrice/Dropbox/college/Astro99/DATA/"

def cogal_downsampled_demo():

    print "loading data: ..."
    cogal = getdata(data_path+'COGAL_all_mom.fits', memmap=True)

    cogal_downsampled = cogal[::4, ::4, ::4]

    print "computing dendrogram: ..."
    d = astrodendro.Dendrogram.compute(
        cogal_downsampled.transpose(2,0,1), 
        min_value=0.01, min_delta=0.005, min_npix=30, verbose=True)

    d.viewer(galactic=True)

    metadata = {}
    metadata['data_unit'] = u.Jy # This is a LIE but Kelvin don't work yet

    catalog = astrodendro.ppv_catalog(d, metadata)

    return d, catalog

    
