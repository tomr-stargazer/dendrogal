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

def transpose_and_downsample_header(input_header, downsample_factor=4, 
                                    transpose_tuple=(2,0,1)):

    df = downsample_factor
    tt = transpose_tuple

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

    new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df + 1
    new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df + 1
    new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df + 1

    return new_header
    

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

    cogal_downsampled_transposed = cogal[::df, ::df, ::df].transpose(2,0,1)

    # beam solid angle. Beam width is 1/8 degree
    omega_B = np.pi * (1/8 * 0.5 * u.deg)**2
    frequency = 115 * u.GHz
    K_to_Jy = u.K.to(u.Jy, 
                     equivalencies=u.brightness_temperature(omega_B, frequency))
    # Convert the data from kelvin to jansky
    cogal_dt_jansky = cogal_downsampled_transposed * K_to_Jy

    cogal_dt_header = transpose_and_downsample_header(cogal_header, df, (2,0,1))
    cogal_dt_wcs = wcs.wcs.WCS(cogal_dt_header)

    v_scale = cogal_header['cdelt1'] * df
    v_unit = u.km / u.s
    l_scale = cogal_header['cdelt2'] * df
    b_scale = cogal_header['cdelt3'] * df

    assert v_scale == cogal_dt_header['cdelt3']
    assert l_scale == cogal_dt_header['cdelt1']
    assert b_scale == cogal_dt_header['cdelt2']

    print "computing dendrogram: ..."
    d = astrodendro.Dendrogram.compute(
        cogal_dt_jansky, 
        min_value=0.01, min_delta=0.005, min_npix=2000//df**3, verbose=True)

    d.viewer(galactic=True)

    metadata = {}
    metadata['data_unit'] = u.Jy # Now this is actually true because I converted!
    metadata['spatial_scale'] = b_scale * u.deg
    metadata['velocity_scale'] = v_scale * v_unit
    metadata['vaxis'] = 0 # keep it this way if you think the input data is (l, b, v)
    metadata['wcs'] = cogal_dt_wcs

    catalog = astrodendro.ppv_catalog(d, metadata)

    return d, catalog, cogal_dt_header, metadata
