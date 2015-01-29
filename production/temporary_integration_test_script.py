""" A temporary hacked-together way to make sure things are working right. """

from load_and_process_data import load_data, permute_data_to_standard_order
from compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog

datacube, header = load_data("DHT03_RCrA_mom.fits")

# _p means permuted
datacube_p, header_p = permute_data_to_standard_order(datacube, header)

d = compute_dendrogram(datacube_p, header_p, min_value=0.1, min_delta=0.1, min_npix=10)

catalog, metadata = compute_catalog(d, header_p)

