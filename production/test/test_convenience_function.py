"""
Tests for convenience_function.py

should run with py.test

"""

import os
import pickle

from numpy.testing import assert_allclose, assert_equal, assert_array_equal, assert_raises
import numpy as np

from astropy.table import Table
from astropy.io.fits import getdata, getheader, Header, writeto
import astropy.units as u
from astrodendro import Dendrogram

from ..convenience_function import save_dendrogram_catalog_output, filename_generator

filepath = "production/test/"

def test_filename_generator():

    # 1. test that no arguments gives an error
    assert_raises(ValueError, filename_generator, [], {})

    # 2. test with some integers
    expected_1 = {}
    expected_1_base = 'a.fits_1.000_2.000_3'
    expected_1['d'] = filepath+expected_1_base+"_d.hdf5"
    expected_1['catalog'] = filepath+expected_1_base+"_catalog.fits"
    expected_1['header'] = filepath+expected_1_base+"_header.fits"
    expected_1['metadata'] = filepath+expected_1_base+"_metadata.p"

    kwargs_1 = {'data_filename':'a.fits',
              'min_value': 1, 
              'min_delta': 2, 
              'min_npix': 3,
              'savepath': filepath}

    assert_equal(filename_generator(**kwargs_1), expected_1)

    # 3. test with some reasonable floats
    expected_2_base = 'a.fits_0.123_12.778_200'
    expected_2_d = filepath+expected_2_base+"_d.hdf5"

    kwargs_2 = {'data_filename':'a.fits',
              'min_value': 0.123, 
              'min_delta': 12.778, 
              'min_npix': 200,
              'savepath': filepath}

    assert_equal(filename_generator(**kwargs_2)['d'], expected_2_d)

    # 4. test with some possibly borderline floats
    expected_3_base = 'a.fits_0.001_3.142_5'
    expected_3_d = filepath+expected_3_base+"_d.hdf5"

    kwargs_3 = {'data_filename':'a.fits',
              'min_value': 0.0009, 
              'min_delta': np.pi, 
              'min_npix': 5,
              'savepath': filepath}

    assert_equal(filename_generator(**kwargs_3)['d'], expected_3_d)

    # 5. test with some negative values
    expected_4_base = 'a.fits_-1.000_-0.550_-5'
    expected_4_d = filepath+expected_4_base+"_d.hdf5"

    kwargs_4 = {'data_filename':'a.fits',
              'min_value': -1, 
              'min_delta': -0.55, 
              'min_npix': -5,
              'savepath': filepath}

    assert_equal(filename_generator(**kwargs_4)['d'], expected_4_d)

    # 6. some None's
    expected_5_base = 'a.fits_None_None_None'
    expected_5_d = filepath+expected_5_base+"_d.hdf5"

    kwargs_5 = {'data_filename':'a.fits',
              'min_value': None, 
              'min_delta': None, 
              'min_npix': None,
              'savepath': filepath}

    assert_equal(filename_generator(**kwargs_5)['d'], expected_5_d)

def test_save_dendrogram_catalog_output():

    # 1. construct a silly dendrogram, catalog, header, and metadata
    data = np.zeros((3,3,3))
    data[1,1,1] = 1
    d = Dendrogram.compute(data, min_value=0)

    catalog = Table()
    catalog['test_column'] = np.zeros(10)

    header = Header.fromstring("SIMPLE  =                    T / conforms to FITS standard                      BITPIX  =                   64 / array data type                                NAXIS   =                    3 / number of array dimensions                     NAXIS1  =                    6                                                  NAXIS2  =                    5                                                  NAXIS3  =                    4                                                  CTYPE1  = 'VELO-LSR'                                                            CTYPE2  = 'GLON-CAR'                                                            CTYPE3  = 'GLAT-CAR'                                                            CRVAL1  = ''                                                                    CRVAL2  = ''                                                                    CRVAL3  = ''                                                                    CDELT1  = ''                                                                    CDELT2  = ''                                                                    CDELT3  = ''                                                                    CRPIX1  = ''                                                                    CRPIX2  = ''                                                                    CRPIX3  = ''                                                                    END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ")

    metadata = {}
    metadata['data_unit'] = u.K

    # 2. run save_dendrogram_catalog_output on them - to some safe location
    kwargs = {'data_filename': 'not_real_data.fits', 
              'min_value': 1, 
              'min_delta': 2, 
              'min_npix': 3,
              'savepath': filepath}
    save_dendrogram_catalog_output(d, catalog, header, metadata, **kwargs)

    # 3. reload those things 
    filename_base = filepath+"not_real_data.fits_1.000_2.000_3"

    d2 = Dendrogram.load_from(filename_base+"_d.hdf5")
    catalog2 = Table.read(filename_base+"_catalog.fits")
    header2 = getheader(filename_base+"_header.fits")
    metadata2 = pickle.load(open(filename_base+"_metadata.p", 'rb'))

    #4 and assert they equal the mock things.
    assert_equal(d2.index_map, d.index_map)
    assert_array_equal(catalog2, catalog)
    assert_equal(header2, header)
    assert_equal(metadata2, metadata)

    #5 then delete the saved objects.
    os.remove(filename_base+"_d.hdf5")
    os.remove(filename_base+"_catalog.fits")
    os.remove(filename_base+"_header.fits")
    os.remove(filename_base+"_metadata.p")

def test_reload_dendrogram_catalog_output():
    pass

# def reload_dendrogram_catalog_output(**kwargs):

#     filename_dict = filename_generator(**kwargs)

#     d = astrodendro.Dendrogram.load_from(filename_dict['d'])
#     catalog = astropy.table.Table.read(filename_dict['catalog'])
#     header = getheader(filename_dict['header'])
#     metadata = pickle.load(open(filename_dict['metadata'], 'rb'))

#     return d, catalog, header, metadata
