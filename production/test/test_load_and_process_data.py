""" 
Tests for load_and_process_data.py

should run with py.test

"""

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from matplotlib.cbook import flatten

from astropy.io.fits import getdata, getheader, Header, writeto

from ..load_and_process_data import load_data, permute_data_to_standard_order
from ..config import data_path

def create_test_data():

    filename = "production/test/test_data.fits"

    datacube = np.arange(120).reshape(4,5,6)

    key_bases = ['naxis', 'ctype', 'crval', 'cdelt', 'crpix']

    keys = list(flatten([[x+'1', x+'2', x+'3'] for x in key_bases]))

    header = Header.fromkeys(keys)
    header['ctype1'] = 'GLON-CAR'
    header['ctype2'] = 'GLAT-CAR'
    header['ctype3'] = 'VELO-LSR'

    writeto(open(filename, 'wb'), datacube, header=header)

# maybe a little overkill... basically just testing the `getdata` function
def test_load_data():

    expected_datacube = np.arange(120).reshape(4,5,6)

    expected_header = Header.fromstring("SIMPLE  =                    T / conforms to FITS standard                      BITPIX  =                   64 / array data type                                NAXIS   =                    3 / number of array dimensions                     NAXIS1  =                    6                                                  NAXIS2  =                    5                                                  NAXIS3  =                    4                                                  CTYPE1  = 'GLON-CAR'                                                            CTYPE2  = 'GLAT-CAR'                                                            CTYPE3  = 'VELO-LSR'                                                            CRVAL1  = ''                                                                    CRVAL2  = ''                                                                    CRVAL3  = ''                                                                    CDELT1  = ''                                                                    CDELT2  = ''                                                                    CDELT3  = ''                                                                    CRPIX1  = ''                                                                    CRPIX2  = ''                                                                    CRPIX3  = ''                                                                    END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ")

    datacube, header = load_data("test_data.fits", data_path="production/test/")

    assert_allclose(datacube, expected_datacube)
    assert_equal(header, expected_header)


def test_permute_data_to_standard_order():

    expected_datacube = np.arange(120).reshape(4,5,6).transpose(2,0,1)

    expected_header = Header.fromstring("SIMPLE  =                    T / conforms to FITS standard                      BITPIX  =                   64 / array data type                                NAXIS   =                    3 / number of array dimensions                     NAXIS1  =                    5                                                  NAXIS2  =                    4                                                  NAXIS3  =                    6                                                  CTYPE1  = 'GLAT-CAR'                                                            CTYPE2  = 'VELO-LSR'                                                            CTYPE3  = 'GLON-CAR'                                                            CRVAL1  = ''                                                                    CRVAL2  = ''                                                                    CRVAL3  = ''                                                                    CDELT1  = ''                                                                    CDELT2  = ''                                                                    CDELT3  = ''                                                                    CRPIX1  = ''                                                                    CRPIX2  = ''                                                                    CRPIX3  = ''                                                                    END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ")
    
    # if something's wrong with load_data then it should fail the previous test
    datacube, header = permute_data_to_standard_order(*load_data("test_data.fits", data_path="production/test/"))

    assert_allclose(datacube, expected_datacube)
    assert_equal(header, expected_header)




