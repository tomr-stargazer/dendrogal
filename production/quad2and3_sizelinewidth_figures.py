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
