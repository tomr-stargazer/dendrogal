"""
This is a "dummy generator" so I don't have to recompute/reload certain things
over and over when working on my code.

"""

from dendrogal.first_quadrant_cloud_extraction import *

firstquad_output = first_quad_dendrogram()
d, catalog, header, metadata = firstquad_output

cloud_catalog = export_firstquad_catalog(firstquad_output)

