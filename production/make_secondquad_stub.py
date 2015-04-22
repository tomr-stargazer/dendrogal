""" 
Loads and computes the second quadrant dendrogram, so that it can be imported. 

Since the dendrogram parameters are (at some fundamental level) now fixed, 
having this stub allows us to run & vary the cloud extraction many times
without having to reload the dendrogram & catalog itself, allowing us to
avoid a time-bottleneck.

"""

from __future__ import division

from astrodendro_analysis.production.convenience_function import load_permute_dendro_catalog

data_filename = "DHT17_Quad2_bw_mominterp.fits"
dendrogram_kwargs = {'min_value' : 0.31/2,
                     'min_delta' : 0.31,
                     'min_npix' : 20}

d, catalog, header, metadata = load_permute_dendro_catalog(data_filename, **dendrogram_kwargs)
