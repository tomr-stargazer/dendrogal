from __future__ import division

from astropy.io.fits import getdata
from astropy import wcs

import astrodendro


def demo(filename):

    min_delta = 1
    min_value = 1
    min_npix = 1000

    datacube, datacube_header = getdata(filename, memmap=True,
                                        header=True)    
    datacube_wcs = wcs.wcs.WCS(datacube_header)

    datacube_wcs.wcs.bounds_check(pix2world=False)

    d = astrodendro.Dendrogram.compute(
        datacube,
        min_value=min_value, min_delta=min_delta,  #these are arbitrary
        min_npix=min_npix, verbose=True, wcs=datacube_wcs)

    return d