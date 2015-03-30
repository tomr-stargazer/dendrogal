""" 
Combines loading and computing dendrograms into one function.

I am implementing caching options to save even more time. 

"""

import os.path
from functools import wraps
import pickle
import datetime

# import astropy.io.fits as fits
import astropy
from astropy.io.fits import getheader

import astrodendro

from .load_and_process_data import load_data, permute_data_to_standard_order
from .compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog


def memoize(func):
    """ Decorator function, ripped off of astrodendro.analysis.memoize()"""

    @wraps(func)
    def wrapper(filename, **kwargs):        
        beginning = datetime.datetime.now()
        try:
            # LOAD THE dendrogram & catalog FROM WHEREVER IT WOULD LIVE IF IT WERE SAVED
            return reload_dendrogram_catalog_output(data_filename=filename, **kwargs)
        except IOError:
            # generate and save the thing, then return it
            output = func(filename=filename, **kwargs)
            save_dendrogram_catalog_output(*output, data_filename=filename, **kwargs)
            return output
        finally:
            end = datetime.datetime.now()
            time_elapsed = (end - beginning)
            print " *** Dendrogram+catalog loading/generation took {0}".format(time_elapsed)


    return wrapper


@memoize
def load_permute_dendro_catalog(filename, min_value=None, min_delta=None, min_npix=None):

    datacube, header = permute_data_to_standard_order(*load_data(filename))

    d = compute_dendrogram(datacube, header, min_value=min_value, min_delta=min_delta, min_npix=min_npix)
    catalog, metadata = compute_catalog(d, header)

    return d, catalog, header, metadata 

# def load_permute_dendro_catalog_memoized(filename, **kwargs):
#     try:
#         # LOAD THE dendrogram & catalog FROM WHEREVER IT WOULD LIVE IF IT WERE SAVED
#         return reload_dendrogram_catalog_output(data_filename=filename, **kwargs)
#     except IOError:
#         # generate and save the thing, then return it
#         output = load_permute_dendro_catalog(filename=filename, **kwargs)
#         save_dendrogram_catalog_output(*output, data_filename=filename, **kwargs)
#         return output

savepath = os.path.expanduser("~/Documents/Code/astrodendro_analysis/production/saved_dendrograms/")

def save_dendrogram_catalog_output(d, catalog, header, metadata, **kwargs):

    filename_dict = filename_generator(**kwargs)

    d.save_to(filename_dict['d'])
    catalog.write(filename_dict['catalog'], overwrite=True)
    header.tofile(filename_dict['header'], clobber=True)
    pickle.dump(metadata, open(filename_dict['metadata'], 'wb'))

    return None


def reload_dendrogram_catalog_output(**kwargs):

    filename_dict = filename_generator(**kwargs)

    d = astrodendro.Dendrogram.load_from(filename_dict['d'])
    catalog = astropy.table.Table.read(filename_dict['catalog'])
    header = getheader(filename_dict['header'])
    metadata = pickle.load(open(filename_dict['metadata'], 'rb'))

    return d, catalog, header, metadata


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
