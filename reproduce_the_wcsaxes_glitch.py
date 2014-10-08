"""
Snippet of plotting code that reproduces my WCSAxes glitch.

"""

from __future__ import division

import os.path
import matplotlib.pyplot as plt

from astropy import wcs
from astropy.io.fits import getdata

from wcsaxes import WCSAxes

data_path = os.path.expanduser("~/Dropbox/College/Astro99/DATA/") # Alter this if you're not Tom Rice

data_file = 'COGAL_local_mom.fits'


def downsample_and_transpose_data_and_header(input_data, input_header, 
                                             downsample_factor=4, 
                                             transpose_tuple=(2,0,1), resample=False, recenter=True):
    
    df = downsample_factor
    tt = transpose_tuple

    if resample and df > 1:
        new_data = resample_3d(input_data, df).transpose(*tt)
    else:
        # Someday I may improve this with something less crude.
        new_data = input_data[::df, ::df, ::df].transpose(*tt)

    new_header = input_header.copy()

    # let's transpose. and downsample.
    new_header['naxis'+str(tt[0]+1)] = input_header['naxis1'] // df
    new_header['naxis'+str(tt[1]+1)] = input_header['naxis2'] // df
    new_header['naxis'+str(tt[2]+1)] = input_header['naxis3'] // df    

    new_header['ctype'+str(tt[0]+1)] = input_header['ctype1']
    new_header['ctype'+str(tt[1]+1)] = input_header['ctype2']
    new_header['ctype'+str(tt[2]+1)] = input_header['ctype3']

    new_header['crval'+str(tt[0]+1)] = input_header['crval1']
    new_header['crval'+str(tt[1]+1)] = input_header['crval2']
    new_header['crval'+str(tt[2]+1)] = input_header['crval3']

    new_header['cdelt'+str(tt[0]+1)] = input_header['cdelt1'] * df
    new_header['cdelt'+str(tt[1]+1)] = input_header['cdelt2'] * df
    new_header['cdelt'+str(tt[2]+1)] = input_header['cdelt3'] * df

    if resample:
        new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df
        new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df
        new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df
    else:
        new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df + 1
        new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df + 1
        new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df + 1

    if recenter:
        return_header = recenter_wcs_header(new_header)
    else:
        return_header = new_header

    return new_data, return_header

def recenter_wcs_header(input_header, central_value=0):
    """ 
    Sets the header CRVAL on zero if it's not already. 

    A quick-and-dirty implementation.

    """

    new_header = input_header.copy()

    # in principle we could iterate through dimensions but.... not yet
    if input_header['CRVAL1'] != central_value:
        # how do we figure out where it hits zero? Use CDELT1 and do some math?
        degrees_correction = central_value - input_header['CRVAL1']
        pixels_correction = degrees_correction / input_header['CDELT1']

        new_header['CRVAL1'] = central_value
        new_header['CRPIX1'] = input_header['CRPIX1'] + pixels_correction

    return new_header


datacube, datacube_header = getdata(data_path+data_file, memmap=True,
	header=True)

datacube_dt, datacube_dt_header = \
	downsample_and_transpose_data_and_header(datacube, datacube_header, 1, (2,0,1), resample=False, recenter=True)
datacube_dt_wcs = wcs.wcs.WCS(datacube_dt_header)

datacube_dt_wcs.wcs.bounds_check(pix2world=False, world2pix=False)

fig = plt.figure(figsize=(13,3.5))
ax_image_limits = [0.1, 0.1, 0.8, 0.8]
ax_image = fig.add_axes(WCSAxes(fig, ax_image_limits, wcs=datacube_dt_wcs, slices=('x','y', 0)))

image = ax_image.imshow(datacube_dt[50, :,:], origin='lower', interpolation='nearest', cmap=plt.cm.gray)

fig.show()