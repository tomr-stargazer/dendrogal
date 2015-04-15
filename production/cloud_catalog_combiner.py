""" This script combines the four quadrants. """

from __future__ import division

import numpy as np
import astropy

from astrodendro_analysis.production.new_cloud_extractor_q1 import export_firstquad_catalog
from astrodendro_analysis.production.second_quadrant_cloud_extraction import export_secondquad_catalog
from astrodendro_analysis.production.third_quadrant_cloud_extraction import export_thirdquad_catalog
from astrodendro_analysis.production.new_cloud_extractor_q4 import export_fourthquad_catalog
from astrodendro_analysis.production.new_cloud_extractor_carina import export_carina_catalog

def extract_and_combine_catalogs(args=[None]*4):

    # The keyword exists so that if you have already computed
    # a dendrogram and catalog for any of the quadrants,
    # you can pass its `d, catalog, header, metadata` tuple
    # in as an element of `args` and skip the longish computation.

    first_cat = export_firstquad_catalog()
    second_cat = export_secondquad_catalog(args[1])
    third_cat = export_thirdquad_catalog(args[2])
    fourth_cat = export_fourthquad_catalog()
    carina_cat = export_carina_catalog()

    # in the combined catalog, we might wanna consider skipping the idx column, or at least renaming it

    first_cat['quadrant'] = 1*np.ones(len(first_cat))
    second_cat['quadrant'] = 2*np.ones(len(second_cat))
    third_cat['quadrant'] = 3*np.ones(len(third_cat))
    fourth_cat['quadrant'] = 4*np.ones(len(fourth_cat))
    carina_cat['quadrant'] = 4*np.ones(len(carina_cat))

    total_table = astropy.table.vstack([first_cat, second_cat, third_cat, fourth_cat, carina_cat])

    return total_table