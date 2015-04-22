""" 
multipanel_catalog_measurements.py

Makes a multi-panel plot. Not meant to be particularly interactive / dynamic.

Front-loads the data as needed. 

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

from astrodendro_analysis.production.catalog_measurement import size_linewidth_slope, cumulative_massfunction_fit, truncated_cloudmass_function
from astrodendro_analysis.production.cloud_catalog_combiner import extract_and_combine_catalogs

all_catalog = extract_and_combine_catalogs()

def extract_inner_galaxy(catalog):

    inner_galaxy = catalog['R_gal'] < 8

    inner_galaxy_catalog = catalog[inner_galaxy]

    return inner_galaxy_catalog


def extract_outer_galaxy(catalog):

    outer_galaxy = catalog['R_gal'] > 9

    outer_galaxy_catalog = catalog[outer_galaxy]

    return outer_galaxy_catalog


def multipanel_size_linewidth():

    fig = plt.figure(figsize=(8,8))

    overall = fig.add_subplot(331)
    plt.setp(overall.get_xticklabels(), visible=False)
    overall.text(1, 1.3, "All clouds", fontsize=18)

    size_array = all_catalog['size']
    linewidth_array = all_catalog['v_rms']

    overall.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    outer_N = fig.add_subplot(337, sharex=overall, sharey=overall)
    outer_N.set_xlabel(r"$\log_{10}(R/\rm{pc})$")
    outer_N.text(1, 1.3, "Outer (North)", fontsize=18)    

    N = all_catalog['x_cen'] < 180
    N_catalog = all_catalog[N]

    size_array = extract_outer_galaxy(N_catalog)['size']
    linewidth_array = extract_outer_galaxy(N_catalog)['v_rms']

    outer_N.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    inner_N = fig.add_subplot(334, sharex=overall, sharey=overall)
    plt.setp(inner_N.get_xticklabels(), visible=False)
    inner_N.text(1, 1.3, "Inner (North)", fontsize=18)    

    size_array = extract_inner_galaxy(N_catalog)['size']
    linewidth_array = extract_inner_galaxy(N_catalog)['v_rms']

    inner_N.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')


    outer_S = fig.add_subplot(338, sharex=overall, sharey=overall)
    outer_S.set_xlabel(r"$\log_{10}(R/\rm{pc})$")    
    outer_S.text(1, 1.3, "Outer (South)", fontsize=18)        
    plt.setp(outer_S.get_yticklabels(), visible=False)

    S = all_catalog['x_cen'] > 180
    S_catalog = all_catalog[S]

    size_array = extract_outer_galaxy(S_catalog)['size']
    linewidth_array = extract_outer_galaxy(S_catalog)['v_rms']

    outer_S.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    inner_S = fig.add_subplot(335, sharex=overall, sharey=overall)
    inner_S.text(1, 1.3, "Inner (South)", fontsize=18)        
    plt.setp(inner_S.get_xticklabels(), visible=False)
    plt.setp(inner_S.get_yticklabels(), visible=False)

    size_array = extract_inner_galaxy(S_catalog)['size']
    linewidth_array = extract_inner_galaxy(S_catalog)['v_rms']

    inner_S.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    outer_combined = fig.add_subplot(339, sharex=overall, sharey=overall)
    plt.setp(outer_combined.get_yticklabels(), visible=False)
    outer_combined.set_xlabel(r"$\log_{10}(R/\rm{pc})$")
    outer_combined.text(1, 1.3, "Outer (all)", fontsize=18)            

    size_array = extract_outer_galaxy(all_catalog)['size']
    linewidth_array = extract_outer_galaxy(all_catalog)['v_rms']

    outer_combined.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    inner_combined = fig.add_subplot(336, sharex=overall, sharey=overall)
    inner_combined.text(1, 1.3, "Inner (all)", fontsize=18)                
    plt.setp(inner_combined.get_xticklabels(), visible=False)
    plt.setp(inner_combined.get_yticklabels(), visible=False)

    size_array = extract_inner_galaxy(all_catalog)['size']
    linewidth_array = extract_inner_galaxy(all_catalog)['v_rms']

    inner_combined.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    overall.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5])
    overall.set_xticks(np.linspace(0.5,2.5, 25), minor=True)
    overall.set_xlim([0.7, 2.6])

    overall.set_yticks(np.linspace(-0, 1.5, 4))
    overall.set_yticks(np.linspace(-0, 1.5, 4*5), minor=True)
    overall.set_ylim([-0.3, 1.6])

    plt.tight_layout(pad=0.5)

    return fig

