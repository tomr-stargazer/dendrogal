"""
Plots the output of run_noise_experiment.multiple_noise_trials_experiment.

Assumes the output has the following structure:

the multiple noise trials experiment gives an output dict

In [12]: output.keys()
Out[12]:
['0.27:3',
 '0.27:2',
 ...
 '0.18:2']

In [14]: output['0.09:0']
Out[14]:
(0,
 0.09,
 408,
 94077069.217720449,
 {'inner_larson_A': 0.50936507099353556,
  'inner_larson_beta': 0.50532497596646464,
  'outer_larson_A': 0.036268556012489767,
  'outer_larson_beta': 0.97668782390932773},
 {'inner_M0': 9290523.4,
  'inner_N0': 3.4717549,
  'inner_gamma': -1.717338,
  'outer_M0': 3894510.8,
  'outer_N0': 0.0,
  'outer_gamma': -1.7741693})


in other words: for each trial, return a dict of 
noise_added: trial number

containing a tuple

(i, noise_level, n_clouds, mass, extract_result['larson'], extract_result['mspec'])


 """

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy.table


import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def build_table_from_output_dict(output_dict):
    """ Traverses the output dict to build a plottable table. """

    table = astropy.table.Table()

    column_names = (
        "trial_number",
        "noise_added",
        "n_clouds",
        "total_mass",
        'inner_larson_A',
        'inner_larson_beta',
        'outer_larson_A',
        'outer_larson_beta',
        'inner_M0',
        'inner_N0',
        'inner_gamma',
        'outer_M0',
        'outer_N0',
        'outer_gamma')

    for name in column_names:

        col = astropy.table.Column(data=[], name=name)

        table.add_column(col)

    for key, value in output_dict.items():

        # we are gonna try and add things row by row

        preamble = value[0:4]
        larson, mspec = value[-2:]

        larson_names = ['inner_larson_A', 'inner_larson_beta',
                        'outer_larson_A', 'outer_larson_beta']
        larson_tuple = tuple(larson[name] for name in larson_names)

        mspec_names = ['inner_M0', 'inner_N0', 'inner_gamma',
                       'outer_M0', 'outer_N0', 'outer_gamma']
        mspec_tuple = tuple(mspec[name] for name in mspec_names)

        table.add_row( (preamble + larson_tuple + mspec_tuple))

    return table


def plot_noise_experiment(output_table):
    """
    Makes a six-panel plot of the output of the noise trials.

    """

    # "real" values from Quadrant I for comparison - these are plotted as dotted lines
    real_nclouds = 389
    real_totalmass = 9.69

    real_larson_A = 0.52
    err_larson_A = 0.08
    real_larson_beta = 0.51
    err_larson_beta = 0.04

    real_mspec_M0 = 8.19
    err_mspec_M0 = 2.57
    real_mspec_gamma = -1.59
    err_mspec_gamma = 0.13

    if type(output_table) is dict:
        output_table = build_table_from_output_dict(output_table)

    fig = plt.figure(figsize=(8,8))

    ax_nclouds = fig.add_subplot(321)
    ax_totalmass = fig.add_subplot(322, sharex=ax_nclouds)

    ax_larson_A = fig.add_subplot(323, sharex=ax_nclouds)
    ax_larson_beta = fig.add_subplot(324, sharex=ax_nclouds)

    ax_mspec_M0 = fig.add_subplot(325, sharex=ax_nclouds)
    ax_mspec_gamma = fig.add_subplot(326, sharex=ax_nclouds)

    noise_added = output_table['noise_added']

    ax_nclouds.plot(noise_added, output_table['n_clouds'], 'ko', ms=5)
    ax_nclouds.set_ylim(0, 450)
    ax_nclouds.set_yticks(np.linspace(0,400,5))
    # ax_nclouds.set_ylabel("Number of clouds")

    ax_nclouds.text(0.025, 200, "(a) Number of clouds", fontsize=14, family='serif')

    ax_nclouds.plot([0, 0.5], [real_nclouds]*2, 'b--', lw=0.5, scalex=False, scaley=False)


    ax_nclouds.set_xticks([0]+list(set(noise_added) - set([0.045])))
    ax_nclouds.set_xlim(0, 0.4)


    ax_totalmass.plot(noise_added, output_table['total_mass']/1e7, 'ko', ms=5)
    ax_totalmass.set_ylim(0, 11)
    # ax_totalmass.set_ylabel("Total mass / $10^7 M_\odot$")

    ax_totalmass.text(0.025, 3, "(b) Total mass of clouds\n\n$(\\times 10^7 M_\odot)$", fontsize=14, family='serif')

    ax_totalmass.plot([0, 0.5], [real_totalmass]*2, 'b--', lw=0.5, scalex=False, scaley=False)


    ax_larson_A.plot(noise_added, output_table['inner_larson_A'], 'ko', ms=5)
    ax_larson_A.set_ylim(0, 0.65)

    ax_larson_A.plot([0, 0.5], [real_larson_A]*2, 'b--', lw=0.5, scalex=False, scaley=False)
    ax_larson_A.plot([0, 0.5], [real_larson_A+err_larson_A]*2, 'k:', lw=0.5, scalex=False, scaley=False)
    ax_larson_A.plot([0, 0.5], [real_larson_A-err_larson_A]*2, 'k:', lw=0.5, scalex=False, scaley=False)

    ax_larson_A.text(0.025, 0.1, "(c) Size-linewidth $A$\nfrom $\\sigma_v  = A \\times R^\\beta$", fontsize=14, family='serif')


    ax_larson_beta.plot(noise_added, output_table['inner_larson_beta'], 'ko', ms=5)
    ax_larson_beta.set_ylim(0, 0.6)

    ax_larson_beta.plot([0, 0.5], [real_larson_beta]*2, 'b--', lw=0.5, scalex=False, scaley=False)
    ax_larson_beta.plot([0, 0.5], [real_larson_beta+err_larson_beta]*2, 'k:', lw=0.5, scalex=False, scaley=False)
    ax_larson_beta.plot([0, 0.5], [real_larson_beta-err_larson_beta]*2, 'k:', lw=0.5, scalex=False, scaley=False)

    ax_larson_beta.text(0.025, 0.1, "(d) Size-linewidth $\\beta$\nfrom $\\sigma_v  = A \\times R^\\beta$", fontsize=14, family='serif')


    ax_mspec_M0.plot(noise_added, output_table['inner_M0']/1e6, 'ko', ms=5)
    ax_mspec_M0.set_ylim(0, 12)

    ax_mspec_M0.plot([0, 0.5], [real_mspec_M0]*2, 'b--', lw=0.5, scalex=False, scaley=False)
    ax_mspec_M0.plot([0, 0.5], [real_mspec_M0+err_mspec_M0]*2, 'k:', lw=0.5, scalex=False, scaley=False)
    ax_mspec_M0.plot([0, 0.5], [real_mspec_M0-err_mspec_M0]*2, 'k:', lw=0.5, scalex=False, scaley=False)

    ax_mspec_M0.set_xlabel("Noise added (K)", family='serif', fontsize=14)
    ax_mspec_M0.text(0.025, 1, "(e) Mass spectrum\n\ntruncation mass $M_0$", fontsize=14, family='serif')


    ax_mspec_gamma.plot(noise_added, output_table['inner_gamma'], 'ko', ms=5)
    ax_mspec_gamma.set_ylim(-2, -1.4)

    ax_mspec_gamma.plot([0, 0.5], [real_mspec_gamma]*2, 'b--', lw=0.5, scalex=False, scaley=False)
    ax_mspec_gamma.plot([0, 0.5], [real_mspec_gamma+err_mspec_gamma]*2, 'k:', lw=0.5, scalex=False, scaley=False)
    ax_mspec_gamma.plot([0, 0.5], [real_mspec_gamma-err_mspec_gamma]*2, 'k:', lw=0.5, scalex=False, scaley=False)

    ax_mspec_gamma.text(0.025, -1.9, "(f) Mass spectrum slope $\\gamma$", fontsize=14, family='serif')


    return fig

