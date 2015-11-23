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

from dendrogal.production.convenience_function import load_permute_dendro_catalog, load_permute_data
from dendrogal.production.load_and_process_data import load_data, permute_data_to_standard_order
from dendrogal.dame_moment_masking import moment_mask, gsmooth_cube
from dendrogal.production.compute_dendrogram_and_catalog import compute_dendrogram, compute_catalog
from dendrogal.production.cloud_extractor_q1 import first_quad_cloud_catalog, compile_firstquad_catalog

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

# so for the "part 2" thing above... we'll have to take the INTERPOLATED data and moment-mask it manually yeah?

def compute_noise_added_processed_dendrogram(noise_added=0, original_noise=0.18, smoothed_rms_noise=None):

    if noise_added > 0:

        # get some data
        print "loading data..."
        datacube, header = load_permute_data(data_filename, memmap=False)

        print "HEADER:"
        print header
        print datacube.shape

        # add some noise to it
        print "adding noise to data..."
        noise_cube = np.random.normal(scale=noise_added, size=datacube.shape)
        noise_added_cube = datacube + noise_cube

        # moment-mask it
        print "moment-masking data..."
        total_noise = (noise_added**2+original_noise**2)**(1/2) # unused, I think

        if smoothed_rms_noise is None:
            smoothed_rms_noise = 0.05

        moment_masked_cube = moment_mask(noise_added_cube, total_noise, velocity_smoothing=4, smoothed_rms_noise=smoothed_rms_noise)

        # now we'll dendrogram it
        print "dendrogramming data..."
        d = compute_dendrogram(moment_masked_cube, header, min_value=min_value, min_delta=min_delta, min_npix=min_npix)
        catalog, metadata = compute_catalog(d, header)

        return d, catalog

    else:

        from dendrogal.production.make_firstquad_stub import d, catalog

        return d, catalog

def extract_properties_from_dendrogram_catalog(dendrogram, raw_catalog):

    # we want to extract the following:

    # (Size-linewidth results, mass spectrum results, total mass, number of clouds.)

    processed_catalog = first_quad_cloud_catalog(d=dendrogram, catalog=raw_catalog)
    reduced_catalog = compile_firstquad_catalog(processed_catalog, dendrogram=dendrogram)

    n_clouds = len(reduced_catalog)
    total_mass = np.nansum(reduced_catalog['mass'])

    print "number of clouds: {0}\ntotal mass: {1:.2e} Msun".format(n_clouds, total_mass)

    output_dict = {}

    output_dict['n_clouds'] = n_clouds
    output_dict['total_mass'] = total_mass

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

    raw_rmss = []

    print "We have added {0:.5f} noise.".format(noise_added)
    print ""
    print "Here are some stats on the raw cube:"

    for (x, y) in coordinate_pairs:

        raw_spectrum = datacube[:, x, y]
        rms = np.nanstd(raw_spectrum)

        print "  rms at {0}, {1}: {2:.4f}".format(x, y, rms)

        raw_rmss.append(rms)

    print "Mean raw rms: {0:.4f}".format(np.mean(raw_rmss))


    smooth_rmss = []

    print "Here are some stats on the smooth cube:"

    for (x, y) in coordinate_pairs:

        smooth_spectrum = smooth_cube[:, x, y]
        rms = np.nanstd(smooth_spectrum)

        print "  rms at {0}, {1}: {2:.4f}".format(x, y, rms)

        smooth_rmss.append(rms)

    print "Mean smooth rms: {0:.4f}".format(np.mean(smooth_rmss))

    return np.mean(smooth_rmss)


def multiple_noise_trials_experiment():

    # ok so here's the story:

    # first get the expected noise to use for moment-masking

    # then do the whole shebang and get your d, catalog

    # then 

    noise_levels = [0.05, 0.1, 0.2]

    n_times_per_noise_level = 3

    output_dict = {}

    for noise_level in noise_levels:

        smoothed_rms_noise = 1.2*inspect_noise_in_first_quadrant_smoothed_cube_before_masking(noise_added=noise_level)

        for i in range(n_times_per_noise_level):

            d, catalog = compute_noise_added_processed_dendrogram(noise_added=noise_level, smoothed_rms_noise=smoothed_rms_noise)

            extract_result = extract_properties_from_dendrogram_catalog(d, catalog)
            n_clouds, mass = extract_result['n_clouds'], extract_result['total_mass']

            output_dict['{0}:{1}'.format(noise_level, i)] = (n_clouds, mass)

    return output_dict
