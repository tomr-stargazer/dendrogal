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
from dendrogal.production.plot_catalog_measurements import plot_size_linewidth_with_nearfar, plot_cmf_except_farU

all_catalog = extract_and_combine_catalogs()


def multipanel_size_linewidth():

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

            plt.text(0.9, 1.1, sky_text_dict[sky_region]+gal_text_dict[galactic_region], fontsize=14)

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

            if subplot_number == 9:
                ax.legend(bbox_to_anchor=(0.45, 0.3), bbox_transform=fig.transFigure)



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


def multipanel_cmf():
    """
    Makes a multi-panel cloud mass function plot.

    """

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

            ax = fig.add_subplot(3, 3, subplot_number, 
                sharex=share_ax, sharey=share_ax)

            fig.axes_list.append(ax)

            catalog = all_catalog[sky_subset_dict[sky_region] & gal_subset_dict[galactic_region]]

            if galactic_region == 'outer':
                max_mass = 1e6
                min_mass = 1e4
            else:
                max_mass = 1e7
                min_mass = 1e5

            cmf_output = plot_cmf_except_farU(catalog, ax, labels=False, max_mass=max_mass, min_mass=min_mass, bins=30, hist_range=(3.7,7))

            M_0, N_0, gamma = cmf_output[0]

            M_0_power = int(np.floor(np.log10(M_0)))
            M_0_coeff = M_0 / (10**M_0_power)

            plt.text(4.5, 3e2, sky_text_dict[sky_region]+gal_text_dict[galactic_region], fontsize=14)

            fit_string = ("$\gamma = {{{0:.2f}}}$\n"
                          "$M_0 = {{{1:.1f}}} \\times 10^{{{2}}}$".format(gamma, M_0_coeff, M_0_power))

            plt.text(3.9, 0.75, fit_string, fontsize=14)

            if sky_to_column_dict[sky_region] != 1:
                plt.setp(ax.get_yticklabels(), visible=False)
            else:
                ax.set_ylabel("n(M > M')")

            if (gal_to_row_dict[galactic_region] == 3):
                pass
            elif (gal_to_row_dict[galactic_region] == 2):
                if (sky_to_column_dict[sky_region] == 3):
                    plt.setp(ax.get_xticklabels(), visible=False)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            if subplot_number in ( 9,):
                ax.set_xlabel(r"log (M$_{GMC}$ / M$_\odot$)")

            # if subplot_number == 9:
            #     ax.legend(bbox_to_anchor=(0.45, 0.3), bbox_transform=fig.transFigure)


    # there's a bug with tight_layout in mac os x, hence all the draw()s.
    # workaround from https://github.com/matplotlib/matplotlib/issues/2654
    fig.canvas.draw()
    fig.tight_layout(pad=0.5)
    fig.canvas.draw() 

    return fig
