"""
This is a quick module to ensure that compute_flux treats Kelvin units properly.

"""

from __future__ import division

import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs
from astropy.io.fits import getdata

from demo import downsample_and_transpose_data_and_header

data_path = "/Users/tsrice/Dropbox/college/Astro99/DATA/"

def compare_jansky_to_kelvin(data_file, downsample_factor=4, transpose_tuple=(2,0,1),
                             min_value=0.01, min_delta=0.005, min_npix=1000):

    df = downsample_factor
    tt = transpose_tuple

    print "loading data: ..."
    datacube, datacube_header = getdata(data_path+data_file, memmap=True,
                                        header=True)

    print "transposing, downsampling, and unit-converting data: ..."
    datacube_dt, datacube_dt_header = \
      downsample_and_transpose_data_and_header(datacube, datacube_header, df, tt)
    datacube_dt_wcs = wcs.wcs.WCS(datacube_dt_header)

    beam_size = 1/8 * u.deg

    # Convert the data from kelvin to jansky-per-beam
    omega_beam = np.pi * (0.5 * beam_size)**2 # Beam width is 1/8 degree
    frequency = 115 * u.GHz
    K_to_Jy = u.K.to(u.Jy, equivalencies=
                     u.brightness_temperature(omega_beam, frequency))
    datacube_dt_jansky_perbeam = datacube_dt * K_to_Jy 

    print "computing dendrogram: ..."
    d_jansky = astrodendro.Dendrogram.compute(
    	datacube_dt_jansky_perbeam,
        min_value=min_value*K_to_Jy, min_delta=min_delta*K_to_Jy,  #these are arbitrary
        min_npix=min_npix//df**3, verbose=True)

    d_kelvin = astrodendro.Dendrogram.compute(
        datacube_dt,
        min_value=min_value, min_delta=min_delta,  #these are arbitrary
        min_npix=min_npix//df**3, verbose=True)

    print len(list(d_jansky.all_structures))
    print len(list(d_kelvin.all_structures))

    v_scale = datacube_dt_header['cdelt3']
    v_unit = u.km / u.s
    l_scale = datacube_dt_header['cdelt1']
    b_scale = datacube_dt_header['cdelt2']
    
    metadata_jansky = {}
    metadata_jansky['data_unit'] = u.Jy / u.beam # According to A. Ginsburg
    metadata_jansky['spatial_scale'] = b_scale * u.deg
    metadata_jansky['velocity_scale'] = v_scale * v_unit
    metadata_jansky['wavelength'] = (c.c / frequency).to('mm')
    metadata_jansky['beam_major'] = beam_size
    metadata_jansky['beam_minor'] = beam_size    
    metadata_jansky['vaxis'] = 0 # keep it this way if you think the (post-downsample/transposed) input data is (l, b, v)
    metadata_jansky['wcs'] = datacube_dt_wcs

    metadata_kelvin = metadata_jansky.copy()
    metadata_kelvin['data_unit'] = u.K # overwrite Jy/beam

    catalog_jansky = astrodendro.ppv_catalog(d_jansky, metadata_jansky)
    catalog_kelvin = astrodendro.ppv_catalog(d_kelvin, metadata_kelvin)

    # if catalog['flux'].unit.is_equivalent('Jy'):
    #     # Workaround because flux is computed wrong

    #     flux = quantify_column(catalog['flux'])
    #     area_exact = catalog['area_exact'].unit*catalog['area_exact'].data

    #     flux_kelvin = flux.to('K', equivalencies=u.brightness_temperature(area_exact, frequency))

    #     flux_kelvin_kms_deg2 = flux_kelvin * metadata['velocity_scale'] * area_exact

    #     catalog.add_column(astropy.table.Column(data=flux_kelvin_kms_deg2, name='flux_kelvin_kms_deg2'))

    return_dict = {}
    return_dict['jansky'] = (d_jansky, catalog_jansky)
    return_dict['kelvin'] = (d_kelvin, catalog_kelvin)

    return return_dict
