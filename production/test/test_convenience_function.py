"""
Tests for convenience_function.py

should run with py.test

"""

import os
import pickle

from numpy.testing import assert_allclose, assert_equal, assert_array_equal
import numpy as np

from astropy.table import Table
from astropy.io.fits import getdata, getheader, Header, writeto
import astropy.units as u
from astrodendro import Dendrogram

from ..convenience_function import save_dendrogram_catalog_output

filepath = "production/test/"

def test_filename_generator():
    pass

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
