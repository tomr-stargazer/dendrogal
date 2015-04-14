"""
Distance disambiguator.

Calculates the "probability" of a given cloud's distance being correct,
given its latitude and its size and linewidth.

The basic approach is to calculate a two-tailed p-value, like this:
https://en.wikipedia.org/wiki/One-_and_two-tailed_tests

and I explicitly assume normal distributions.

"""

from __future__ import division

import pdb

import numpy as np
from scipy.stats import norm

import astropy.units as u

from astrodendro_analysis.production.calculate_distance_dependent_properties import assign_properties


def p_given_sigmas(N_sigmas):

    p = 1 - (norm.cdf(N_sigmas) - norm.cdf(-N_sigmas))

    return p


def p_from_size_linewidth(R, sigma_v, A_coefficient, B_coefficient, scatter_in_log10_R=0.265):
    """
    Probability that the chosen R is consistent with the sigma_v given A & B.

    Assumes a power law of the form:

    sigma_v = A * R**B

    and that the scatter of log10(R) is Gaussian around 1/B*log10(sigma_v/A)

    """

    log_R_offset = np.log10(R) - (1 / B_coefficient) * np.log10(sigma_v / A_coefficient)

    N_sigmas = np.abs(log_R_offset / scatter_in_log10_R)

    return p_given_sigmas(N_sigmas)


def p_from_latitude(z, molecular_HWHM_height=60*u.pc):
    """
    Probability that the chosen z is consistent with the molecular scale height.

    Assumes that molecular clouds follow a Gaussian distribution in z.

    Conversion from half-width-half-maximum to sigmas comes from here:
    http://mathworld.wolfram.com/GaussianFunction.html

    """

    # sqrt(2 ln 2) is about 1.1774
    sigma = molecular_HWHM_height / np.sqrt(2*np.log(2))

    N_sigmas = np.abs(z / sigma)

    return p_given_sigmas(N_sigmas)


def calculate_p_nearfar(catalog, return_intermediates=False,
                        A_coefficient=0.269502604526,
                        B_coefficient=0.602835495068,
                        scatter_in_log10_R=0.265,
                        molecular_HWHM_height=60*u.pc
                        ):

    near_catalog = catalog.copy(copy_data=True)
    near_catalog['distance'] = catalog['near_distance']
    far_catalog = catalog.copy(copy_data=True)
    far_catalog['distance'] = catalog['far_distance']

    assign_properties(near_catalog)
    assign_properties(far_catalog)

    p_near_larson = p_from_size_linewidth(near_catalog['size'], near_catalog['v_rms'],
                                          A_coefficient, B_coefficient)
    p_far_larson = p_from_size_linewidth(far_catalog['size'], far_catalog['v_rms'],
                                         A_coefficient, B_coefficient)

    p_near_latitude = p_from_latitude((near_catalog['z_gal']).to(u.pc),
                                      molecular_HWHM_height=molecular_HWHM_height)

    p_far_latitude = p_from_latitude((far_catalog['z_gal']).to(u.pc),
                                     molecular_HWHM_height=molecular_HWHM_height)

    p_near = p_near_larson * p_near_latitude
    p_far = p_far_larson * p_far_latitude

    if return_intermediates:
        near_output = [near_catalog, p_near_larson, p_near_latitude, p_near]
        far_output = [far_catalog, p_far_larson, p_far_latitude, p_far]
        return near_output, far_output
    else:
        return p_near, p_far


def distance_disambiguator(catalog, **kwargs):
    """
    Chooses the best distance based on size-linewidth fit & latitude.

    """

    p_near, p_far = calculate_p_nearfar(catalog, False, **kwargs)

    best_distance = np.zeros_like(catalog['near_distance'])

    use_near_distance = (p_near >= p_far)
    use_far_distance = (p_far > p_near)

    best_distance[use_near_distance] = catalog['near_distance'][use_near_distance]
    best_distance[use_far_distance] = catalog['far_distance'][use_far_distance]

    pdb.set_trace()

    return best_distance
