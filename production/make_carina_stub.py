""" 
Loads and computes the Carina dendrogram, so that it can be imported. 

Since the dendrogram parameters are (at some fundamental level) now fixed, 
having this stub allows us to run & vary the cloud extraction many times
without having to reload the dendrogram & catalog itself, allowing us to
avoid a time-bottleneck.

"""

from __future__ import division

from dendrogal.production.convenience_function import load_permute_dendro_catalog

min_value = 0.17
min_delta = 0.17
min_npix = 20

data_filename = "DHT33_Carina_mominterp.fits"
dendrogram_kwargs = {'min_value': min_value,
                     'min_delta': min_delta,
                     'min_npix': min_npix}

d, catalog, header, metadata = load_permute_dendro_catalog(data_filename, **dendrogram_kwargs)
