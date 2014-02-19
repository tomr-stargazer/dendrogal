"""
This is a demo to open, downsample, and dendrogram some sample data.

The data are allowed to be the "actual" thing we're trying to dendrogram I guess.

"""

from __future__ import division

import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs
from astropy.io.fits import getdata

data_path = "/Users/tsrice/Dropbox/college/Astro99/DATA/"

def downsample_and_transpose_data_and_header(input_data, input_header, 
                                             downsample_factor=4, 
                                             transpose_tuple=(2,0,1)):
    
    df = downsample_factor
    tt = transpose_tuple

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

    new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df + 1
    new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df + 1
    new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df + 1

    return new_data, new_header
    

def cogal_downsampled_demo(downsample_factor=4, transpose_tuple=(2,0,1)):

    if type(downsample_factor) != int:
        raise TypeError("Downsample factor must be an Integer")
    elif downsample_factor < 1:
        raise ValueError("Downsample factor must be positive")
    elif downsample_factor > 100:
        raise ValueError("Downsample factor must be less than 100")
    
    df = downsample_factor
    tt = transpose_tuple

    print "loading data: ..."
    cogal, cogal_header = getdata(data_path+'COGAL_all_mom.fits', memmap=True,
                                  header=True)

    print "transposing, downsampling, and unit-converting data: ..."
    cogal_dt, cogal_dt_header = \
      downsample_and_transpose_data_and_header(cogal, cogal_header, df, tt)      
    cogal_dt_wcs = wcs.wcs.WCS(cogal_dt_header)

    beam_size = 1/8 * u.deg

    # Convert the data from kelvin to jansky-per-beam
    omega_beam = np.pi * (0.5 * beam_size)**2 # Beam width is 1/8 degree
    frequency = 115 * u.GHz
    K_to_Jy = u.K.to(u.Jy, equivalencies=
                     u.brightness_temperature(omega_beam, frequency))
    cogal_dt_jansky_perbeam = cogal_dt * K_to_Jy 

    print "computing dendrogram: ..."
    d = astrodendro.Dendrogram.compute(
        cogal_dt_jansky_perbeam, 
        min_value=0.01*K_to_Jy, min_delta=0.005*K_to_Jy, 
        min_npix=2000//df**3, verbose=True)

    d.viewer(galactic=True)

    v_scale = cogal_dt_header['cdelt3']
    v_unit = u.km / u.s
    l_scale = cogal_dt_header['cdelt1']
    b_scale = cogal_dt_header['cdelt2']
    
    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam # According to A. Ginsburg
    metadata['spatial_scale'] = b_scale * u.deg
    metadata['velocity_scale'] = v_scale * v_unit
    metadata['wavelength'] = (c.c / frequency).to('mm')
    metadata['beam_major'] = beam_size
    metadata['beam_minor'] = beam_size    
    metadata['vaxis'] = 0 # keep it this way if you think the input data is (l, b, v)
    metadata['wcs'] = cogal_dt_wcs

    catalog = astrodendro.ppv_catalog(d, metadata)

    return d, catalog, cogal_dt_header, metadata
