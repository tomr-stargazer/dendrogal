""" 
Combines loading and computing dendrograms into one function.

I am implementing caching options to save even more time. 

"""

import os.path
from functools import wraps
import pickle
# import astropy.io.fits as fits
import astropy
from astropy.io.fits import getheader

import astrodendro

from .load_and_process_data import load_data, permute_data_to_standard_order
from .compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog

def load_permute_dendro_catalog(filename, min_value=None, min_delta=None, min_npix=None):

    datacube, header = permute_data_to_standard_order(*load_data(filename))

    d = compute_dendrogram(datacube, header, min_value=min_value, min_delta=min_delta, min_npix=min_npix)
    catalog, metadata = compute_catalog(d, header)

    return d, catalog, header, metadata 

savepath = os.path.expanduser("~/Documents/Code/astrodendro_analysis/production/saved_dendrograms/")

def save_dendrogram_catalog_output(d, catalog, header, metadata, **kwargs):

    filename_dict = filename_generator(**kwargs)

    d.save_to(filename_dict['d'])
    catalog.write(filename_dict['catalog'], overwrite=True)
    header.tofile(filename_dict['header'], clobber=True)
    pickle.dump(metadata, open(filename_dict['metadata'], 'wb'))

    return None


def filename_generator(data_filename=None, min_value=None, min_delta=None, min_npix=None, savepath=savepath):

    filename_dict = {}

    if data_filename is None: 
        raise ValueError("`data_filename` must be provided!")

    if min_value is None:
        min_value_str = str(min_value)
    else:
        min_value_str = '{:.3f}'.format(min_value)

    if min_delta is None:
        min_delta_str = str(min_delta)
    else:
        min_delta_str = '{:.3f}'.format(min_delta)

    min_npix_str = str(min_npix)

    filename_base = "{0}_{1}_{2}_{3}".format(data_filename, min_value_str, min_delta_str, min_npix_str)

    filename_dict['d'] = savepath+filename_base+"_d.hdf5"
    filename_dict['catalog'] = savepath+filename_base+"_catalog.fits"
    filename_dict['header'] = savepath+filename_base+"_header.fits"
    filename_dict['metadata'] = savepath+filename_base+"_metadata.p"

    return filename_dict


