"""
Makes the size-linewidth figures for the second quadrant.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt
import astropy.table

from dendrogal.production.cloud_extractor_q2 import export_secondquad_catalog
from dendrogal.production.cloud_extractor_q2 import d as quad2_d

from dendrogal.production.cloud_extractor_q3 import export_thirdquad_catalog
from dendrogal.production.cloud_extractor_q3 import d as quad3_d

from dendrogal.production.catalog_measurement import size_linewidth_slope
from dendrogal.production.plot_catalog_measurements import plot_size_linewidth_with_nearfar

output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

second_quadrant_catalog = export_secondquad_catalog()
third_quadrant_catalog = export_thirdquad_catalog()

outer_catalog = astropy.table.vstack([second_quadrant_catalog, third_quadrant_catalog])

odr_outer = size_linewidth_slope(outer_catalog)

def outer_galaxy_size_linewidth_with_residuals(): 

    fig = plt.figure()

    second_larson = fig.add_subplot(221)
    third_larson = fig.add_subplot(222)

    second_residuals_plot = fig.add_subplot(223)
    third_residuals_plot = fig.add_subplot(224)

    second_larson.errorbar(
        second_quadrant_catalog['size'], 
        second_quadrant_catalog['v_rms'],
        xerr=(second_quadrant_catalog['error_size_minus'], 
              second_quadrant_catalog['error_size_plus']),
        fmt='o')
    second_larson.set_xscale('log')
    second_larson.set_yscale('log')

    second_fit = size_linewidth_slope(second_quadrant_catalog)
    second_fit_coefficient = second_fit.beta[0]
    second_fit_exponent = second_fit.beta[1]

    second_residuals = second_quadrant_catalog['v_rms'] - second_fit_coefficient * (second_quadrant_catalog['size'])**(second_fit_exponent)

    second_residuals_plot.errorbar(
        second_residuals,        
        second_quadrant_catalog['v_rms'],
        xerr=(second_quadrant_catalog['error_size_minus'], 
              second_quadrant_catalog['error_size_plus']),
        fmt='o'
        )
    second_residuals_plot.set_yscale('log')

    print "second quadrant residual RMS: {0}".format(np.std(second_residuals))


    third_larson.errorbar(
        third_quadrant_catalog['size'], 
        third_quadrant_catalog['v_rms'],
        xerr=(third_quadrant_catalog['error_size_minus'], 
              third_quadrant_catalog['error_size_plus']),
        fmt='o')
    third_larson.set_xscale('log')
    third_larson.set_yscale('log')

    third_fit = size_linewidth_slope(third_quadrant_catalog)
    third_fit_coefficient = third_fit.beta[0]
    third_fit_exponent = third_fit.beta[1]

    third_residuals = third_quadrant_catalog['v_rms'] - third_fit_coefficient * (third_quadrant_catalog['size'])**(third_fit_exponent)

    third_residuals_plot.errorbar(
        third_residuals,        
        third_quadrant_catalog['v_rms'],
        xerr=(third_quadrant_catalog['error_size_minus'], 
              third_quadrant_catalog['error_size_plus']),
        fmt='o'
        )
    third_residuals_plot.set_yscale('log')

    print "third quadrant residual RMS: {0}".format(np.std(third_residuals))


    return fig


def new_multipanel_outer_galaxy(allcat):

    fig = plt.figure(figsize=(7,7))

    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2, sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(2,2,3, sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(2,2,4, sharex=ax1, sharey=ax1)

    ax_list = [ax4, ax2, ax1, ax3]

    special_allcat_1 = allcat[(allcat['v_cen']<-20 ) & (allcat['x_cen'] < 80) & (allcat['x_cen'] > 20) ]
    special_allcat_2 = allcat[(np.abs(allcat['v_cen'])>20 ) & (allcat['x_cen'] < 160) & (allcat['x_cen'] > 80) ]
    special_allcat_3 = allcat[(np.abs(allcat['v_cen'])>20 ) & (allcat['x_cen'] > 200) & (allcat['x_cen'] < 280) ]
    special_allcat_4 = allcat[(allcat['v_cen']>20 ) & (allcat['x_cen'] < 340) & (allcat['x_cen'] > 280) ]

    subcat_list = [special_allcat_1, special_allcat_2, special_allcat_3, special_allcat_4]    

    name_list = ['IQ', 'IIQ', 'IIIQ', 'IVQ']

    for ax, subcat, name in zip(ax_list, subcat_list, name_list):
        size_linewidth_output = plot_size_linewidth_with_nearfar(subcat, ax, labels=False, distance_threshold=20)

        fit_coefficient = size_linewidth_output.beta[0]
        sd_coefficient = size_linewidth_output.sd_beta[0]
        fit_exponent = size_linewidth_output.beta[1]
        sd_exponent = size_linewidth_output.sd_beta[1]

        fit_string = "$A = {{{0:.2f} \pm {2:.2f}}}$,\n$\\beta ={{{1:.2f} \pm {3:.2f}}}$".format(fit_coefficient, fit_exponent, sd_coefficient, sd_exponent) 

        ax.text(1.85, -0.2, fit_string, fontsize=12)
        ax.text(0.9, 1.1, name, fontsize=16, family='serif')

    ax1.set_ylabel(r"$\log(\sigma_v)$", fontsize=16)
    ax3.set_ylabel(r"$\log(\sigma_v)$", fontsize=16)

    ax3.set_xlabel(r"$\log(R/\rm{pc})$", fontsize=16)
    ax4.set_xlabel(r"$\log(R/\rm{pc})$", fontsize=16)

    ax1.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax1.set_xticks(np.linspace(0.5,3, 30), minor=True)

    ax1.set_yticks(np.linspace(-0, 1.5, 7))
    ax1.set_yticks(np.linspace(-0, 1.5, 4*5), minor=True)

    ax1.set_xlim(0.75, 3)
    ax1.set_ylim(-0.25, 1.25)

    return fig


def save_multipanel_outer_galaxy(allcat):
    fig = new_multipanel_outer_galaxy(allcat)
    fig.savefig(output_path+"outer_galaxy_fourpanel_larson.pdf", bbox_inches='tight')


def new_singlepanel_outer_galaxy(allcat):

    fig = plt.figure(figsize=(7,7))

    ax1 = fig.add_subplot(1,1,1)

    quad_1_outer_criteria = (allcat['v_cen']<-20 ) & (allcat['x_cen'] < 80) & (allcat['x_cen'] > 20) 
    quad_2_outer_criteria = (np.abs(allcat['v_cen'])>20 ) & (allcat['x_cen'] < 160) & (allcat['x_cen'] > 80) 
    quad_3_outer_criteria = (np.abs(allcat['v_cen'])>20 ) & (allcat['x_cen'] > 200) & (allcat['x_cen'] < 280) 
    quad_4_outer_criteria = (allcat['v_cen']>20 ) & (allcat['x_cen'] < 340) & (allcat['x_cen'] > 280) 

    outer_galaxy_allcat = allcat[quad_2_outer_criteria | quad_3_outer_criteria]

    size_linewidth_output = plot_size_linewidth_with_nearfar(outer_galaxy_allcat, ax1, labels=False, distance_threshold=20)

    ax1.set_xlim(0.8, 2.5)
    ax1.set_ylim(-0.4, 1)

    fit_coefficient = size_linewidth_output.beta[0]
    sd_coefficient = size_linewidth_output.sd_beta[0]
    fit_exponent = size_linewidth_output.beta[1]
    sd_exponent = size_linewidth_output.sd_beta[1]

    fit_string = "$A = {{{0:.2f} \pm {2:.2f}}}$,\n$\\beta ={{{1:.2f} \pm {3:.2f}}}$".format(fit_coefficient, fit_exponent, sd_coefficient, sd_exponent) 

    ax1.text(1.85, -0.1, fit_string, fontsize=24)
    ax1.text(0.9, 1.1, "IIQ + IIIQ", fontsize=18, family='serif')

    ax1.set_ylabel(r"$\log(\sigma_v)$", fontsize=20)
    ax1.set_xlabel(r"$\log(R/\rm{pc})$", fontsize=20)

    ax1.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax1.set_xticks(np.linspace(0.5,3, 30), minor=True)

    ax1.set_yticks(np.linspace(-0, 1.5, 7))
    ax1.set_yticks(np.linspace(-0, 1.5, 4*5), minor=True)

    ax1.set_xlim(0.75, 3)
    ax1.set_ylim(-0.25, 1.25)

    return fig


def save_singlepanel_outer_galaxy(allcat):
    fig = new_singlepanel_outer_galaxy(allcat)
    fig.savefig(output_path+"outer_galaxy_onepanel_larson.pdf", bbox_inches='tight')


def new_multipanel_inner_galaxy(allcat):

    fig = plt.figure(figsize=(7,7))

    # ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3, sharex=ax2, sharey=ax2)
    ax4 = fig.add_subplot(2,2,4, sharex=ax2, sharey=ax2)

    ax_list = [ax4, ax3, ax2]

    quad_1_inner_criteria = (allcat['v_cen']>20 ) & (allcat['x_cen'] < 80) & (allcat['x_cen'] > 20)
    quad_4_inner_criteria = (allcat['v_cen']<-20 ) & (allcat['x_cen'] < 340) & (allcat['x_cen'] > 280) 

    special_allcat_1 = allcat[quad_1_inner_criteria]
    special_allcat_4 = allcat[quad_4_inner_criteria]
    special_allcat_1_4 = allcat[quad_1_inner_criteria | quad_4_inner_criteria]

    subcat_list = [special_allcat_1, special_allcat_4, special_allcat_1_4]    

    name_list = ['IQ', 'IVQ', 'IQ + IVQ']

    for ax, subcat, name in zip(ax_list, subcat_list, name_list):
        size_linewidth_output = plot_size_linewidth_with_nearfar(
            subcat, ax, labels=False, distance_threshold=30, 
            alternate_style=True)

        fit_coefficient = size_linewidth_output.beta[0]
        sd_coefficient = size_linewidth_output.sd_beta[0]
        fit_exponent = size_linewidth_output.beta[1]
        sd_exponent = size_linewidth_output.sd_beta[1]

        fit_string = "$A = {{{0:.2f} \pm {2:.2f}}}$,\n$\\beta ={{{1:.2f} \pm {3:.2f}}}$".format(fit_coefficient, fit_exponent, sd_coefficient, sd_exponent) 

        ax.text(1.85, -0.2, fit_string, fontsize=12)
        ax.text(0.9, 1.1, name, fontsize=16, family='serif')


    ax3.legend(bbox_to_anchor=(0.49, 0.65), bbox_transform=fig.transFigure, numpoints=1)

    ax3.set_ylabel(r"$\log(\sigma_v)$", fontsize=16)
    ax2.set_ylabel(r"$\log(\sigma_v)$", fontsize=16)

    ax4.set_xlabel(r"$\log(R/\rm{pc})$", fontsize=16)
    ax3.set_xlabel(r"$\log(R/\rm{pc})$", fontsize=16)

    ax3.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax3.set_xticks(np.linspace(0.5,3, 30), minor=True)

    ax3.set_yticks(np.linspace(-0, 1.5, 7))
    ax3.set_yticks(np.linspace(-0, 1.5, 4*5), minor=True)

    ax3.set_xlim(0.75, 3)
    ax3.set_ylim(-0.25, 1.25)

    fig.ax3 = ax3
    fig.ax4 = ax4

    return fig


def save_multipanel_inner_galaxy(allcat):
    fig = new_multipanel_inner_galaxy(allcat)
    fig.savefig(output_path+"inner_galaxy_threepanel_larson.pdf", bbox_inches='tight')


