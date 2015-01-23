""" Saves us some time. """

from .load_and_process_data import load_data, permute_data_to_standard_order
from .compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog

def load_permute_dendro_catalog(filename, min_value=None, min_delta=None, min_npix=None):

    datacube, header = permute_data_to_standard_order(*load_data(filename))

    d = compute_dendrogram(datacube, header, min_value=min_value, min_delta=min_delta, min_npix=min_npix)
    catalog, metadata = compute_catalog(d, header)

    return d, catalog, header, metadata 
