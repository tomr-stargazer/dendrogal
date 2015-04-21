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

    first_cat['survey'] = 8*np.ones(len(first_cat))
    second_cat['survey'] = 17*np.ones(len(second_cat))
    third_cat['survey'] = 31*np.ones(len(third_cat))
    fourth_cat['survey'] = 36*np.ones(len(fourth_cat))
    carina_cat['survey'] = 33*np.ones(len(carina_cat))

    total_table = astropy.table.vstack([first_cat, second_cat, third_cat, fourth_cat, carina_cat])

    return total_table

def print_results_by_quadrant(total_table=None):

    if total_table is None:
        total_table = extract_and_combine_catalogs()

    first_cat = total_table[total_table['quadrant']==1]
    second_cat = total_table[total_table['quadrant']==2]
    third_cat = total_table[total_table['quadrant']==3]
    fourth_cat = total_table[total_table['survey']==36]
    carina_cat = total_table[total_table['survey']==33]
    combined_fourth_cat = total_table[total_table['quadrant']==4]

    print "Quadrant       |  N_clouds  | Mass  "
    print "------------------------------------"
    print "I               | {0:4d}     | {1:.2e} ".format(len(first_cat), np.sum(first_cat['mass']))
    print "II              | {0:4d}     | {1:.2e} ".format(len(second_cat), np.sum(second_cat['mass']))
    print "III             | {0:4d}     | {1:.2e} ".format(len(third_cat), np.sum(third_cat['mass']))
    print "IV (no Carina)  | {0:4d}     | {1:.2e} ".format(len(fourth_cat), np.sum(fourth_cat['mass']))
    print "IV (Carina-only)| {0:4d}     | {1:.2e} ".format(len(carina_cat), np.sum(carina_cat['mass']))
    print "IV (combined)   | {0:4d}     | {1:.2e} ".format(len(combined_fourth_cat), np.sum(combined_fourth_cat['mass']))
    print "All combined    | {0:4d}     | {1:.2e} ".format(len(total_table), np.sum(total_table['mass']))

    return total_table
