"""
Figuring out why I get `nan` in the properties of structure 1517 & friends.

See https://github.com/dendrograms/astrodendro/issues/122

"""

from __future__ import division

import pickle
import os.path
import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c
from astropy import wcs

from astrodendro.analysis import PPVStatistic, ScalarStatistic

header = astropy.io.fits.header.Header.fromfile("debug_structure_1517_header.fits")
header_wcs = wcs.wcs.WCS(header)

v_scale = header['cdelt3']
v_unit = u.km / u.s
l_scale = header['cdelt1']
b_scale = header['cdelt2']

beam_size = 1/8 * u.deg
frequency = 115 * u.GHz

metadata = {}
metadata['data_unit'] = u.K
metadata['spatial_scale'] = b_scale * u.deg
metadata['velocity_scale'] = v_scale * v_unit
metadata['wavelength'] = frequency # formerly: (c.c / frequency).to('mm') but now compute_flux can handle frequency in spectral equivalency
metadata['beam_major'] = beam_size
metadata['beam_minor'] = beam_size    
metadata['vaxis'] = 0 # keep it this way if you think the (post-downsample/transposed) input data is (l, b, v)
metadata['wcs'] = header_wcs

#In [67]: metadata
#Out[67]:
# {'beam_major': <Quantity 0.125 deg>,
#  'beam_minor': <Quantity 0.125 deg>,
#  'data_unit': Unit("K"),
#  'spatial_scale': <Quantity 0.25 deg>,
#  'vaxis': 0,
#  'velocity_scale': <Quantity 2.6 km / s>,
#  'wavelength': <Quantity 115.0 GHz>,
#  'wcs': <astropy.wcs.wcs.WCS at 0x10673a638>}

#In [63]: d.data.shape
#Out[63]: 
shape_tuple = (246, 40, 1440)
metadata['shape_tuple'] = shape_tuple

def properties_1517():
    # In [63]: d, catalog, header, metadata = cogal_deep_resampled_demo(resample=2, min_npix=20)

    # In [61]: d[1517].indices()
    # Out[61]:
    indices = (np.array([124, 124, 124, 123, 123, 123, 123, 123, 123, 123, 123, 124, 124,
                         123, 123, 123, 124, 123, 124, 124, 124, 123, 124, 124, 124]),
               np.array([ 9, 11, 10, 10, 12,  8,  8,  9,  9,  9, 10,  8,  8, 11, 11, 10, 10,
                          9, 10, 10,  9,  8,  9,  9,  9]),
               np.array([   2, 1438, 1439,    0, 1439,    2,    1,    2, 1439,    0,    1,
                            2,    1, 1439,    0, 1439,    1,    1,    0, 1438, 1438,    3,
                            0,    1, 1439]))

    # In [62]: d[1517].values()
    # Out[62]:
    values = np.array([ 0.7245028 ,  0.67389417,  1.1429944 ,  4.07517517,  0.69078934,
                       0.7946701 ,  0.55045485,  1.23698175,  0.60745633,  1.6095072 ,
                       1.14116812,  0.64040899,  0.87998176,  1.24877799,  1.62031388,
                       1.76993275,  1.27678406,  2.46003664,  3.16056681,  0.65349865,
                       0.69018042,  0.64809525,  1.11399925,  2.93979144,  1.02465391])


    return indices, values

def properties_1451():
    # In [8]: d[1451].indices(subtree=True)
    # Out[8]:
    indices = (np.array([122, 122, 122, 121, 121, 123, 121, 121, 122, 122, 122, 123, 122,
        121, 120, 121, 123, 121, 122, 123, 121, 122, 121, 122, 121, 123,
        121, 121, 120, 122, 122, 121, 121, 121, 123, 121, 123, 121, 120,
        122, 121, 122, 121, 122, 122, 121, 121, 120, 120, 122, 123, 120,
        123, 122, 120, 122, 122, 123, 120, 123, 120, 122, 120, 122, 122,
        121, 123, 121, 121, 123, 121, 123, 121, 122, 122, 122, 122, 122,
        120, 125, 125, 125, 125, 123, 123, 123, 124, 124, 124, 124, 124,
        125, 124, 124, 124, 124, 124, 124, 124, 124, 125, 124, 125, 124,
        124, 124, 123, 123, 123, 123, 123, 123, 123, 123, 124, 124, 123,
        123, 123, 124, 123, 124, 124, 124, 123, 124, 124, 124, 124, 125,
        125, 125, 124, 124, 124, 125, 124, 125, 125, 125, 124, 125, 125,
        125, 124, 125, 125, 124, 125, 124, 124, 125]),
    np.array([15, 10,  9, 11, 12, 13, 13, 11, 11, 12, 14, 12, 13,  9, 10, 10, 11,
        12, 12, 12, 11, 15,  9, 12, 11, 15, 10, 14, 12, 15, 14, 12, 12, 10,
        12, 11, 15, 10, 12, 14,  9,  9, 13, 11, 12, 13,  9, 10, 10, 12, 14,
        11, 14, 10, 10, 12, 14, 14,  9, 10, 10, 13, 11, 11, 13, 10, 11, 11,
        10, 13,  9, 12, 14, 13,  9, 12, 13, 13, 11,  9, 11, 11,  8, 11, 13,
        10, 12,  9, 12, 11,  8,  9, 11,  8, 11, 11,  7,  9,  6, 10, 10,  6,
        10,  9, 11, 10, 10, 12,  8,  8,  9,  9,  9, 10,  8,  8, 11, 11, 10,
        10,  9, 10, 10,  9,  8,  9,  9,  9,  8,  9,  9,  8, 10,  9, 10,  9,
         9, 12,  9, 11, 10,  8, 10, 10,  9, 10, 10, 10, 10,  9, 10, 11]),
    np.array([1439, 1438, 1435, 1436, 1434, 1438, 1435, 1437, 1432, 1437, 1437,
        1435, 1435, 1438,    0, 1437, 1435, 1432, 1435, 1434, 1433, 1437,
        1439, 1434, 1434, 1436, 1433, 1433, 1438, 1438, 1436, 1438, 1437,
        1434, 1433, 1435, 1437, 1436, 1437, 1438, 1435, 1437, 1433, 1434,
        1433, 1432, 1436, 1437, 1438, 1438, 1437, 1438, 1436, 1434, 1439,
        1432, 1435, 1435, 1437, 1435, 1434, 1433, 1434, 1435, 1432, 1438,
        1434, 1438, 1439, 1437, 1437, 1436, 1432, 1437, 1436, 1436, 1438,
        1436, 1437,    6,    0,    3,    7, 1438,    0, 1438,    9,    8,
        1439, 1439,    5,    3,    0,    3,    8,    9,    2,    3,    2,
           7,    7,    3,    5,    2, 1438, 1439,    0, 1439,    2,    1,
           2, 1439,    0,    1,    2,    1, 1439,    0, 1439,    1,    1,
           0, 1438, 1438,    3,    0,    1, 1439,    6,    4,    5,    5,
           8,    5,    6,    8,    7,    9,    7,    8,    5,    6,    3,
           2,    6,    8,    6,    3,    4,    4,    4,    9]))

    values = np.array([ 0.32914662,  0.33257139,  0.32770073,  0.41263199,  0.7063905 ,
        0.37838566,  0.36689389,  0.8697077 ,  0.58348382,  0.92686141,
        1.33066511,  0.59862816,  0.37496078,  1.00106204,  0.37830937,
        1.53218639,  0.76392436,  0.39444327,  0.98325396,  0.55273807,
        0.83736396,  0.64215922,  0.54025722,  0.83028638,  2.21719122,
        0.47473216,  0.3831799 ,  0.4070003 ,  1.03447151,  0.63508177,
        0.59786725,  1.12313139,  0.65814078,  2.3585149 ,  0.53241837,
        0.9351567 ,  0.49893308,  0.42244935,  0.37161231,  0.72944951,
        0.52374268,  0.83386314,  0.35281491,  0.75867331,  0.90836823,
        0.57648218,  0.7899518 ,  0.76255465,  0.41826344,  0.38432145,
        1.18203533,  1.08820009,  1.09238589,  0.39467144,  0.46050096,
        0.85174763,  0.38493049,  0.86567438,  0.39208424,  0.38926828,
        1.0356127 ,  0.42503691,  0.70258522,  0.74345267,  0.51103342,
        1.3906343 ,  1.15874791,  2.00661373,  0.9361459 ,  0.685462  ,
        1.03820038,  0.35220587,  0.40030336,  0.92069697,  0.54976976,
        0.78645086,  0.62305725,  0.52214432,  0.61324   ,  0.42655897,
        0.35106432,  0.51141381,  0.50266206,  0.33964884,  0.44040966,
        0.38622403,  0.46605635,  0.41757882,  0.3394208 ,  0.49482346,
        0.46696973,  0.35882699,  1.29961491,  0.37419975,  0.38736582,
        0.48766983,  0.43911576,  0.51324058,  0.42846131,  0.39718306,
        0.34261703,  0.46346903,  0.48363638,  0.7245028 ,  0.67389417,
        1.1429944 ,  4.07517517,  0.69078934,  0.7946701 ,  0.55045485,
        1.23698175,  0.60745633,  1.6095072 ,  1.14116812,  0.64040899,
        0.87998176,  1.24877799,  1.62031388,  1.76993275,  1.27678406,
        2.46003664,  3.16056681,  0.65349865,  0.69018042,  0.64809525,
        1.11399925,  2.93979144,  1.02465391,  0.55121589,  0.94360387,
        0.5162847 ,  0.90387821,  1.02594757,  1.05623698,  0.89170182,
        0.81567442,  1.320315  ,  0.85410666,  1.98705518,  0.80768371,
        1.39375472,  0.75258505,  2.37327898,  0.58622336,  0.81879473,
        1.4843936 ,  0.53964841,  0.96879423,  1.49656999,  1.41848803,
        1.52411938,  1.124349  ])

    return indices, values

indices, values = properties_1517()
in1451, va1451 = properties_1451()

def compute_the_things(indices=indices, values=values, shape_tuple=shape_tuple, metadata=metadata, statistic=PPVStatistic):

    # make variable `indices` point to a copy in this scope so we don't mutate the original `indices`
    indices = np.copy(indices)

    if shape_tuple is not None:
        for index_array, shape in zip(indices, shape_tuple):
            # catch simple cases where a structure wraps around the image boundary
            i2 = np.where(index_array < shape/2, index_array+shape, index_array)
            if i2.ptp() < index_array.ptp():  # more compact with wrapping. Use this
                index_array[:] = i2

    stat = ScalarStatistic(values, indices)

    stat.mom1()
    ppv_stat = statistic(stat, metadata)

    return stat, ppv_stat

    # row = {}
    # for lbl in fields:
    #     row[lbl] = getattr(stat, lbl)

    # row = dict((lbl, getattr(stat, lbl))
    #            for lbl in fields)
    # row.update(_idx=struct.idx)

    # # first row
    # if result is None:
    #     sorted_row_keys = sorted(row.keys())
    #     try:
    #         result = Table(names=sorted_row_keys,
    #                        dtype=[int if x == '_idx' else float for x in sorted_row_keys])
    #     except TypeError:  # dtype was called dtypes in older versions of Astropy
    #         result = Table(names=sorted_row_keys,
    #                        dtypes=[int if x == '_idx' else float for x in sorted_row_keys])
    #     for k, v in row.items():
    #         try:  # Astropy API change
    #             result[k].unit = _unit(v)
    #         except AttributeError:
    #             result[k].units = _unit(v)

    # # astropy.table.Table should in future support setting row items from
    # # quantities, but for now we need to strip off the quantities
    # new_row = {}
    # for x in row:
    #     if row[x] is not None:  # in Astropy 0.3+ we no longer need to exclude None items
    #         if isinstance(row[x], Quantity):
    #             new_row[x] = row[x].value
    #         else:
    #             new_row[x] = row[x]
    # result.add_row(new_row)    

