""" This script combines the four quadrants. """

from __future__ import division

import numpy as np
import astropy

from dendrogal.production.cloud_extractor_q1 import export_firstquad_catalog
from dendrogal.production.cloud_extractor_q2 import export_secondquad_catalog
from dendrogal.production.cloud_extractor_q3 import export_thirdquad_catalog
from dendrogal.production.cloud_extractor_q4 import export_fourthquad_catalog
from dendrogal.production.cloud_extractor_carina import export_carina_catalog

from dendrogal.production.config import catalog_output_path

def extract_and_combine_catalogs():

    first_cat = export_firstquad_catalog()
    second_cat = export_secondquad_catalog()
    third_cat = export_thirdquad_catalog()
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


def save_catalog_to_file():

    catalog = extract_and_combine_catalogs()

    # reprogram some columns

    catalog['x_cen'].name = 'l (deg)'
    catalog['y_cen'].name = 'b (deg)'
    catalog['v_cen'].name = 'v (km/s)'

    catalog['flux'].name = 'summed_intensity'

    # remove a bunch of columns

    columns_to_remove = [
       '_idx',
    # 'area_ellipse',
    # 'area_exact',
    # 'major_sigma',
    # 'minor_sigma',
    # 'position_angle',
    # 'radius',
    # 'v_cen',
    # # 'v_rms',
    # 'x_cen',
    # 'y_cen',
    # 'flux_true',
    'flux_clipped',
    # 'level',
    # 'n_descendants',
    'is_leaf',
    # 'fractional_gain',
    'on_edge',
    'v_split',
    'max_vsplit' ]
    # 'near_distance',
    # 'error_near_distance_plus',
    # 'error_near_distance_minus',
    # 'far_distance',
    # 'error_far_distance_plus',
    # 'error_far_distance_minus',
    # 'distance',
    # 'error_distance_plus',
    # 'error_distance_minus',
    # 'KDA_resolution',
    # 'p_near',
    # 'p_far',
    # 'size',
    # 'error_size_plus',
    # 'error_size_minus',
    # 'mass',
    # 'error_mass_plus',
    # 'error_mass_minus',
    # 'virial_alpha',
    # 'error_virial_alpha_plus',
    # 'error_virial_alpha_minus',
    # 'pressure',
    # 'error_pressure_plus',
    # 'error_pressure_minus',
    # 'x_sol',
    # 'y_sol',
    # 'z_sol',
    # 'x_gal',
    # 'y_gal',
    # 'z_gal',
    # 'R_gal',
    # 'phi_gal',
    # 'near_size',
    # 'far_size',
    # 'near_z_gal',
    # 'far_z_gal',
    # 'near_mass',
    # 'far_mass',
    # 'quadrant',
    # 'survey']

    for col in columns_to_remove:
        del catalog[col]

    # add an index column
    index_column = astropy.table.Column(data=1+np.arange(len(catalog)), name='cloud index')
    catalog.add_column(index_column, index=0)

    catalog.write(catalog_output_path+"cloud_catalog.fits", format='fits', overwrite=True)


def print_results_by_quadrant(total_table=None):

    if total_table is None:
        total_table = extract_and_combine_catalogs()

    first_cat = total_table[total_table['quadrant']==1]
    second_cat = total_table[total_table['quadrant']==2]
    third_cat = total_table[total_table['quadrant']==3]
    fourth_cat = total_table[total_table['survey']==36]
    carina_cat = total_table[total_table['survey']==33]
    combined_fourth_cat = total_table[total_table['quadrant']==4]

    print "Quadrant        |  N_clouds  | Mass  "
    print "------------------------------------"
    print "I               | {0:4d}     | {1:.2e} ".format(len(first_cat), np.nansum(first_cat['mass']))
    print "II              | {0:4d}     | {1:.2e} ".format(len(second_cat), np.nansum(second_cat['mass']))
    print "III             | {0:4d}     | {1:.2e} ".format(len(third_cat), np.nansum(third_cat['mass']))
    print "IV (no Carina)  | {0:4d}     | {1:.2e} ".format(len(fourth_cat), np.nansum(fourth_cat['mass']))
    print "IV (Carina-only)| {0:4d}     | {1:.2e} ".format(len(carina_cat), np.nansum(carina_cat['mass']))
    print "IV (combined)   | {0:4d}     | {1:.2e} ".format(len(combined_fourth_cat), np.nansum(combined_fourth_cat['mass']))
    print "All combined    | {0:4d}     | {1:.2e} ".format(len(total_table), np.nansum(total_table['mass']))


    inner_cat = total_table[total_table['R_gal'] <= 8.34]
    outer_cat = total_table[total_table['R_gal'] > 8.34]

    inner_1 = first_cat[first_cat['R_gal'] <= 8.34]
    inner_4 = fourth_cat[fourth_cat['R_gal'] <= 8.34]

    print ""
    print "Region          |  N_clouds  | Mass  "
    print "inner           | {0:4d}     | {1:.2e} ".format(len(inner_cat), np.nansum(inner_cat['mass']))
    print "outer           | {0:4d}     | {1:.2e} ".format(len(outer_cat), np.nansum(outer_cat['mass']))
    print "inner I         | {0:4d}     | {1:.2e} ".format(len(inner_1), np.nansum(inner_1['mass']))
    print "inner IV        | {0:4d}     | {1:.2e} ".format(len(inner_4), np.nansum(inner_4['mass']))

    return total_table


def print_ambiguous_results_by_quadrant(total_table):

    first_cat_ambig = total_table[(total_table['quadrant']==1) & (total_table['KDA_resolution']=='A')]
    # no ambiguities in 2/3
    # second_cat_ambig = total_table[(total_table['quadrant']==2) & (total_table['KDA_resolution']=='A')]
    # third_cat_ambig = total_table[(total_table['quadrant']==3) & (total_table['KDA_resolution']=='A')]
    fourth_cat_ambig = total_table[(total_table['survey']==36) & (total_table['KDA_resolution']=='A')]
    # carina_cat_ambig = total_table[(total_table['survey']==33) & (total_table['KDA_resolution']=='A')]
    combined_fourth_cat_ambig = total_table[(total_table['quadrant']==4) & (total_table['KDA_resolution']=='A')]

    all_ambig = total_table[total_table['KDA_resolution'] == 'A']

    print "Results for ambiguous clouds:"
    print "Quadrant       |  N_clouds  | Mass lower limit  "
    print "------------------------------------"
    print "I               | {0:4d}     | {1:.2e} ".format(len(first_cat_ambig), np.sum(first_cat_ambig['near_mass']))
    # print "II              | {0:4d}     | {1:.2e} ".format(len(second_cat_ambig), np.sum(second_cat_ambig['near_mass']))
    # print "III             | {0:4d}     | {1:.2e} ".format(len(third_cat_ambig), np.sum(third_cat_ambig['near_mass']))
    print "IV (no Carina)  | {0:4d}     | {1:.2e} ".format(len(fourth_cat_ambig), np.sum(fourth_cat_ambig['near_mass']))
    # print "IV (Carina-only)| {0:4d}     | {1:.2e} ".format(len(carina_cat_ambig), np.sum(carina_cat_ambig['near_mass']))
    print "IV (combined)   | {0:4d}     | {1:.2e} ".format(len(combined_fourth_cat_ambig), np.sum(combined_fourth_cat_ambig['near_mass']))
    print "All combined    | {0:4d}     | {1:.2e} ".format(len(all_ambig), np.sum(all_ambig['near_mass']))

    print "Results for ambiguous clouds:"
    print "Quadrant       |  N_clouds  | Mass upper limit  "
    print "------------------------------------"
    print "I               | {0:4d}     | {1:.2e} ".format(len(first_cat_ambig), np.sum(first_cat_ambig['far_mass']))
    # print "II              | {0:4d}     | {1:.2e} ".format(len(second_cat_ambig), np.sum(second_cat_ambig['far_mass']))
    # print "III             | {0:4d}     | {1:.2e} ".format(len(third_cat_ambig), np.sum(third_cat_ambig['far_mass']))
    print "IV (no Carina)  | {0:4d}     | {1:.2e} ".format(len(fourth_cat_ambig), np.sum(fourth_cat_ambig['far_mass']))
    # print "IV (Carina-only)| {0:4d}     | {1:.2e} ".format(len(carina_cat_ambig), np.sum(carina_cat_ambig['far_mass']))
    print "IV (combined)   | {0:4d}     | {1:.2e} ".format(len(combined_fourth_cat_ambig), np.sum(combined_fourth_cat_ambig['far_mass']))
    print "All combined    | {0:4d}     | {1:.2e} ".format(len(all_ambig), np.sum(all_ambig['far_mass']))

