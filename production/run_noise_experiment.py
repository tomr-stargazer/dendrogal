"""
This is a module that aims to address how our results vary when random noise is added.

It focuses on the first quadrant, as a representative subset of the data.

Basic ideas:

1. run the analysis un-molested and measure the output. 
(Size-linewidth results, mass spectrum results, total mass, number of clouds.)

2. Create a lightly salted version of the data and repeat the above, just once though.

3. Once that's all established, do it like 20-100 times (maybe just ten?) 
   with 0.5x, 1.0x, and 2.0x the noise.


"""

from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt

from dendrogal.production.convenience_function import load_permute_dendro_catalog, load_permute_data
from dendrogal.production.load_and_process_data import load_data, permute_data_to_standard_order
from dendrogal.dame_moment_masking import moment_mask, gsmooth_cube
from dendrogal.production.compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog
from dendrogal.production.cloud_extractor_q1 import first_quad_cloud_catalog, compile_firstquad_catalog

from dendrogal.production.catalog_measurement import size_linewidth_slope
from dendrogal.production.mspecfit_wrapper import get_mspec_fit

from dendrogal.production.config import data_path

# these are hard-coded because reasons.
min_value = 0.18
min_delta = 0.18
min_npix = 20

dendrogram_kwargs = {'min_value': min_value,
                     'min_delta': min_delta,
                     'min_npix': min_npix}

# btw the noise here is 0.18 K / px / channel
data_filename = "DHT08_Quad1_interp.fits"
raw_filename = "DHT08_Quad1_raw.fits"

# so for the "part 2" thing above... we'll have to take the INTERPOLATED data and moment-mask it manually yeah?

def compute_noise_added_processed_dendrogram(noise_added=0, original_noise=0.18, smoothed_rms_noise=None):

    if noise_added > 0:

        # get some data
        print "loading data..."
        datacube, header = load_permute_data(data_filename, memmap=False)

        # add some noise to it
        print "adding noise to data..."
        noise_cube = np.random.normal(scale=noise_added, size=datacube.shape)
        noise_added_cube = datacube + noise_cube

        # moment-mask it
        print "moment-masking data..."
        total_noise = (noise_added**2+original_noise**2)**(1/2) # unused, I think

        if smoothed_rms_noise is None:
            smoothed_rms_noise = 0.05

        moment_masked_cube = moment_mask(noise_added_cube, total_noise, velocity_smoothing=3, smoothed_rms_noise=smoothed_rms_noise)

        # now we'll dendrogram it
        print "dendrogramming data..."
        d = compute_dendrogram(moment_masked_cube, header, min_value=min_value, min_delta=min_delta, min_npix=min_npix)
        catalog, metadata = compute_catalog(d, header)

        return d, catalog

    else:

        from dendrogal.production.make_firstquad_stub import d, catalog

        return d, catalog

def extract_properties_from_dendrogram_catalog(dendrogram, raw_catalog, i=0):

    # we want to extract the following:

    # (Size-linewidth results, mass spectrum results, total mass, number of clouds.)

    processed_catalog = first_quad_cloud_catalog(d=dendrogram, catalog=raw_catalog)
    reduced_catalog = compile_firstquad_catalog(processed_catalog, dendrogram=dendrogram)

    pruned_catalog = reduced_catalog[
        # on the sky...
        (
        # northernly
        ((reduced_catalog['x_cen'] > 20) & (reduced_catalog['x_cen'] < 160)) |
        # southernly
        ((reduced_catalog['x_cen'] > 200) & (reduced_catalog['x_cen'] < 340)) 
        ) & 
        # in velocity space
        (np.abs(reduced_catalog['v_cen']) > 20)
        ]

    inner_catalog = pruned_catalog[pruned_catalog['R_gal'] < 8]
    outer_catalog = pruned_catalog[pruned_catalog['R_gal'] > 9]

    outer_max_mass = 3e7
    outer_min_mass = 1e4
    inner_max_mass = 3e7
    inner_min_mass = 1e5

    vetted_inner_catalog = inner_catalog[(~np.isnan(inner_catalog['mass'])) &
                                    (inner_catalog['mass']>inner_min_mass) & (inner_catalog['mass']<inner_max_mass)]

    vetted_outer_catalog = outer_catalog[(~np.isnan(outer_catalog['mass'])) &
                                    (outer_catalog['mass']>outer_min_mass) & (outer_catalog['mass']<outer_max_mass)]

    n_clouds = len(reduced_catalog)
    total_mass = np.nansum(reduced_catalog['mass'])

    print "number of clouds: {0}\ntotal mass: {1:.2e} Msun".format(n_clouds, total_mass)

    output_dict = {}

    output_dict['n_clouds'] = n_clouds
    output_dict['total_mass'] = total_mass

    # Size Linewidth results

    inner_larson_output = size_linewidth_slope(inner_catalog)
    outer_larson_output = size_linewidth_slope(outer_catalog)

    larson_dict = {}

    larson_dict['inner_larson_beta'] = inner_larson_output.beta[1]
    larson_dict['inner_larson_A'] = inner_larson_output.beta[0]

    larson_dict['outer_larson_beta'] = outer_larson_output.beta[1]
    larson_dict['outer_larson_A'] = outer_larson_output.beta[0]

    output_dict['larson'] = larson_dict

    # feed stuff to Mass function

    # inner : truncated
    inner_mspec = get_mspec_fit(vetted_inner_catalog, 'inner_', notrunc=0)

    # outer : non-truncated
    outer_mspec = get_mspec_fit(vetted_outer_catalog, 'outer_', notrunc=1)

    mspec_dict = {}

    mspec_dict['inner_N0'] = inner_mspec['N_0']
    mspec_dict['inner_M0'] = inner_mspec['M_0']
    mspec_dict['inner_gamma'] = inner_mspec['gamma']

    mspec_dict['outer_N0'] = outer_mspec['N_0']
    mspec_dict['outer_M0'] = outer_mspec['M_0']
    mspec_dict['outer_gamma'] = outer_mspec['gamma']

    output_dict['mspec'] = mspec_dict

    return output_dict


def inspect_noise_in_first_quadrant_smoothed_cube_before_masking(noise_added=0):

    print "loading data..."
    datacube, header = load_permute_data(data_filename, memmap=False)

    if noise_added > 0:

        datacube += np.random.normal(scale=noise_added, size=datacube.shape)

    velocity_smoothing=4
    spatial_smoothing=2

    smooth_cube = gsmooth_cube(datacube, [velocity_smoothing/2.3548, spatial_smoothing/2.3548, spatial_smoothing/2.3548], kernelsize_mult=4)

    coordinate_pairs = [ (45, 145), (46, 126), (44, 152), (90, 90), (95, 10) ]

    # for (x, y) in coordinate_pairs:

    #     raw_spectrum = datacube[50:-50, x, y]
    #     rms = np.nanstd(raw_spectrum)

    #     print "  rms at {0}, {1}: {2:.4f}".format(x, y, rms)

    #     raw_rmss.append(rms)

    smooth_rmss = []

    print "Here are some stats on the smooth cube:"

    for (x, y) in coordinate_pairs:

        smooth_spectrum = smooth_cube[50:-50, x, y]
        rms = np.nanstd(smooth_spectrum)

        smooth_rmss.append(rms)

    print "Mean smooth rms: {0:.4f}".format(np.mean(smooth_rmss))

    return np.mean(smooth_rmss)


def multiple_noise_trials_experiment():

    # first get the expected noise to use for moment-masking
    # then do the whole shebang and get your d, catalog

    noise_levels = [0.045, 0.09, 0.18, 0.27, 0.36]

    n_times_per_noise_level = 5

    output_dict = {}

    for i, noise_level in enumerate(noise_levels):

        smoothed_rms_noise = inspect_noise_in_first_quadrant_smoothed_cube_before_masking(noise_added=noise_level)

        for i in range(n_times_per_noise_level):

            d, catalog = compute_noise_added_processed_dendrogram(
                noise_added=noise_level, smoothed_rms_noise=smoothed_rms_noise)

            extract_result = extract_properties_from_dendrogram_catalog(d, catalog, i)
            n_clouds, mass = extract_result['n_clouds'], extract_result['total_mass']

            mass_spectrum_dict = {}

            output_dict['{0}:{1}'.format(noise_level, i)] = (i, noise_level, n_clouds, mass, extract_result['larson'], extract_result['mspec'])

    return output_dict


def compare_first_quadrant_moment_masking(smoothed_rms_noise):

    print "loading data..."
    datacube, header = load_permute_data(data_filename, memmap=False)

    noise = 0.18

    moment_masked_cube = moment_mask(datacube, noise, velocity_smoothing=4, smoothed_rms_noise=None)

    actual_moment_masked_cube, masked_header = load_permute_data("DHT08_Quad1_mominterp.fits", memmap=False)

    return moment_masked_cube, actual_moment_masked_cube


reference_cube, reference_header = load_permute_data("DHT08_Quad1_mom.fits", memmap=False)

raw_cube, raw_header = load_permute_data(raw_filename, memmap=False)


def moment_masking_comparison_diagnostic(smoothed_rms_noise=0.05, reference_cube=reference_cube, data_cube=raw_cube):

    moment_cube = moment_mask(raw_cube, rms_noise=0.18, smoothed_rms_noise=smoothed_rms_noise, velocity_smoothing=3, spatial_smoothing=2)

    diff_cube = np.abs(reference_cube - moment_cube)

    fig = plt.figure()

    ax_reference = fig.add_subplot(131)
    ax_moment = fig.add_subplot(132)
    ax_diff = fig.add_subplot(133)

    reference_lv = np.nansum(reference_cube, axis=1)
    reference_lv[(reference_lv < 0) | np.isinf(reference_lv) | np.isnan(reference_lv)] = 0

    moment_lv = np.nansum(moment_cube, axis=1)
    moment_lv[(moment_lv < 0) | np.isinf(moment_lv) | np.isnan(moment_lv)] = 0

    diff_lv = np.nansum(diff_cube, axis=1)
    diff_lv[(diff_lv < 0) | np.isinf(diff_lv) | np.isnan(diff_lv)] = 0

    reference_image = ax_reference.imshow(np.log10(reference_lv+1), origin='lower', 
            interpolation='nearest', cmap='gray', aspect=2.5)

    moment_image = ax_moment.imshow(np.log10(moment_lv+1), origin='lower', 
            interpolation='nearest', cmap='gray', aspect=2.5)

    diff_image = ax_diff.imshow(np.log10(diff_lv+1), origin='lower', 
            interpolation='nearest', cmap='gray', aspect=2.5)

    print ""
    print smoothed_rms_noise
    print "Max diff: {0:.4f}".format(np.nanmax(diff_cube))
    print np.nanmedian(diff_cube)
    print np.nanstd(diff_cube)
    print "{0:.2e}".format(np.nansum(diff_cube))

    ax_reference.set_title("Smoothed RMS noise assumed: {0}".format(smoothed_rms_noise))
    ax_diff.set_xlabel("Diff std: {0:.3f} Diff sum: {1:.2e}".format(np.nanstd(diff_cube), np.nansum(diff_cube)))

    fig.canvas.draw()

    return fig



