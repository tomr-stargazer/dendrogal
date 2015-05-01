""" 
multipanel_catalog_measurements.py

Makes a multi-panel plot. Not meant to be particularly interactive / dynamic.

Front-loads the data as needed. 

"""

from __future__ import division

import pdb

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

from dendrogal.production.catalog_measurement import size_linewidth_slope, cumulative_massfunction_fit, truncated_cloudmass_function
from dendrogal.production.cloud_catalog_combiner import extract_and_combine_catalogs
from dendrogal.production.plot_catalog_measurements import plot_size_linewidth_with_nearfar

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
    # currently coded like a freshman. it's an abomination. computers do more for less.

    fig = plt.figure(figsize=(8,8))

    overall = fig.add_subplot(331)
    plt.setp(overall.get_xticklabels(), visible=False)
    overall.set_ylabel(r"$\log_{10}(\sigma_v)$")    
    overall.text(1, 1.3, "All clouds", fontsize=18)

    size_array = all_catalog['size']
    linewidth_array = all_catalog['v_rms']

    overall.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    pdb.set_trace()
    size_linewidth_output = size_linewidth_slope(all_catalog)
    fit_coefficient = size_linewidth_output.beta[0]
    fit_exponent = size_linewidth_output.beta[1]
    fit_xs = np.logspace(0, 4, 20)
    fit_ys = fit_coefficient * fit_xs ** fit_exponent    
    plt.plot(np.log10(fit_xs), np.log10(fit_ys), 'g--', zorder=0, scalex=False, scaley=False)

    outer_N = fig.add_subplot(337, sharex=overall, sharey=overall)
    outer_N.set_xlabel(r"$\log_{10}(R/\rm{pc})$")
    outer_N.set_ylabel(r"$\log_{10}(\sigma_v)$")
    outer_N.text(1, 1.3, "Outer (North)", fontsize=18)    

    N = all_catalog['x_cen'] < 180
    N_catalog = all_catalog[N]

    size_array = extract_outer_galaxy(N_catalog)['size']
    linewidth_array = extract_outer_galaxy(N_catalog)['v_rms']

    outer_N.plot(np.log10(size_array), np.log10(linewidth_array), 'k.')

    inner_N = fig.add_subplot(334, sharex=overall, sharey=overall)
    plt.setp(inner_N.get_xticklabels(), visible=False)
    inner_N.set_ylabel(r"$\log_{10}(\sigma_v)$")
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

    # there's a bug with tight_layout in mac os x, hence all the draw()s.
    # workaround from https://github.com/matplotlib/matplotlib/issues/2654
    fig.canvas.draw()
    fig.tight_layout(pad=0.5)
    fig.canvas.draw() 

    return fig


def multipanel_size_linewidth2():
    # currently coded like a freshman. it's an abomination. computers do more for less.

    fig = plt.figure(figsize=(8,8))

    sky_to_column_dict = {'north': 1, 'south': 2, 'north_plus_south': 3}
    gal_to_row_dict = {'inner': 1, 'outer': 2, 'inner_plus_outer': 3}

    sky_text_dict = {'north': 'Northern', 'south': 'Southern', 'north_plus_south' : 'All'}
    gal_text_dict = {'inner': ' inner Galaxy', 'outer': ' outer Galaxy', 'inner_plus_outer' : ' Galaxy (combined)'}

    sky_subset_dict = {}
    sky_subset_dict['north'] = all_catalog['x_cen'] < 180
    sky_subset_dict['south'] = all_catalog['x_cen'] > 180
    sky_subset_dict['north_plus_south'] = (sky_subset_dict['north'] | sky_subset_dict['south'])

    gal_subset_dict = {}
    gal_subset_dict['inner'] = all_catalog['R_gal'] < 8
    gal_subset_dict['outer'] = all_catalog['R_gal'] > 9
    gal_subset_dict['inner_plus_outer'] = np.ones(len(all_catalog), dtype='bool')

    fig.axes_list = []

    # the core loop.
    for i, sky_region in enumerate(['north', 'south', 'north_plus_south']):

        for j, galactic_region in enumerate(['inner', 'outer', 'inner_plus_outer']):

            if (galactic_region == 'inner_plus_outer') and (sky_region != 'north_plus_south'):
                continue

            if i+j == 0:
                share_ax = None
            else:
                share_ax = fig.axes_list[0]

            subplot_number = 3*(gal_to_row_dict[galactic_region]-1) + sky_to_column_dict[sky_region]

            ax = fig.add_subplot(3, 3, subplot_number
                , sharex=share_ax, sharey=share_ax)

            fig.axes_list.append(ax)

            catalog = all_catalog[sky_subset_dict[sky_region] & gal_subset_dict[galactic_region]]

            size_linewidth_output = plot_size_linewidth_with_nearfar(catalog, ax, labels=False)

            fit_coefficient = size_linewidth_output.beta[0]
            fit_exponent = size_linewidth_output.beta[1]

            fit_string = "$A = {{{0:.2f}}}$,\n$\\beta ={{{1:.2f}}}$".format(fit_coefficient, fit_exponent)

            plt.text(1,1, sky_text_dict[sky_region]+gal_text_dict[galactic_region], fontsize=14)

            plt.text(2, -0.2, fit_string, fontsize=14)

            if sky_to_column_dict[sky_region] != 1:
                plt.setp(ax.get_yticklabels(), visible=False)
            else:
                ax.set_ylabel(r"$\log(\sigma_v)$", fontsize=16)

            if (gal_to_row_dict[galactic_region] == 3):
                pass
            elif (gal_to_row_dict[galactic_region] == 2):
                if (sky_to_column_dict[sky_region] == 3):
                    plt.setp(ax.get_xticklabels(), visible=False)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            if subplot_number in ( 9,):
                ax.set_xlabel(r"$\log(R/\rm{pc})$", fontsize=16)



    fig.axes_list[0].set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    fig.axes_list[0].set_xticks(np.linspace(0.5,3, 30), minor=True)
    # fig.axes_list[0].set_xlim([0.7, 2.6])

    fig.axes_list[0].set_yticks(np.linspace(-0, 1.5, 4))
    fig.axes_list[0].set_yticks(np.linspace(-0, 1.5, 4*5), minor=True)
    # fig.axes_list[0].set_ylim([-0.3, 1.6])

    fig.axes_list[0].set_xlim(0.75, 3)
    fig.axes_list[0].set_ylim(-0.25, 1.25)

    # there's a bug with tight_layout in mac os x, hence all the draw()s.
    # workaround from https://github.com/matplotlib/matplotlib/issues/2654
    fig.canvas.draw()
    fig.tight_layout(pad=0.5)
    fig.canvas.draw() 

    return fig
