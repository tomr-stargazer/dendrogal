"""
Distance disambiguator.

Calculates the "probability" of a given cloud's distance being correct,
given its latitude and its size and linewidth.

The basic approach is to calculate a two-tailed p-value, like this:
https://en.wikipedia.org/wiki/One-_and_two-tailed_tests

and I explicitly assume normal distributions.

"""

from __future__ import division

import numpy as np
from scipy.stats import norm

import astropy.units as u

from dendrogal.production.calculate_distance_dependent_properties import assign_properties


def p_given_sigmas(N_sigmas):
    """
    Chance to draw a value N_sigmas from the mean of a normal distribution.

    """

    p = 2 * norm.cdf(-np.abs(N_sigmas))

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
                        A_coefficient=0.5,
                        B_coefficient=0.5,
                        scatter_in_log10_R=0.5,
                        molecular_HWHM_height=60*u.pc
                        ):

    near_catalog = catalog.copy(copy_data=True)
    near_catalog['distance'] = catalog['near_distance']
    far_catalog = catalog.copy(copy_data=True)
    far_catalog['distance'] = catalog['far_distance']

    assign_properties(near_catalog)
    assign_properties(far_catalog)

    p_near_larson = p_from_size_linewidth(near_catalog['size'], near_catalog['v_rms'],
                                          A_coefficient, B_coefficient, scatter_in_log10_R=scatter_in_log10_R)
    p_far_larson = p_from_size_linewidth(far_catalog['size'], far_catalog['v_rms'],
                                         A_coefficient, B_coefficient, scatter_in_log10_R=scatter_in_log10_R)

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


def distance_disambiguator(catalog, ambiguous_threshold=0.1, **kwargs):
    """
    Chooses the best distance based on size-linewidth fit & latitude.

    If neither distance is preferred at the 10\% level then don't pick either.

    """

    p_near, p_far = calculate_p_nearfar(catalog, False, **kwargs)

    best_distance = np.zeros_like(catalog['near_distance'])
    error_best_distance_plus = np.zeros_like(catalog['error_near_distance_plus'])
    error_best_distance_minus = np.zeros_like(catalog['error_near_distance_minus'])
    KDA_resolution = np.array(['-']*len(best_distance))

    # if neither distance has a p >~ 10%, then don't pick either.

    same_distances = catalog['near_distance'] == catalog['far_distance']
    ambiguous_distance = (p_near <= ambiguous_threshold) & (p_far <= ambiguous_threshold) & (p_near/p_far < 100) & (p_far/p_near < 100) & ~same_distances
    use_near_distance = (p_near >= p_far) & ~ambiguous_distance
    use_far_distance = (p_far > p_near) & ~ambiguous_distance

    best_distance[use_near_distance] = catalog['near_distance'][use_near_distance]
    best_distance[use_far_distance] = catalog['far_distance'][use_far_distance]
    best_distance[ambiguous_distance] = np.nan
    error_best_distance_plus[use_near_distance] = catalog['error_near_distance_plus'][use_near_distance]
    error_best_distance_plus[use_far_distance] = catalog['error_far_distance_plus'][use_far_distance]
    error_best_distance_plus[ambiguous_distance] = np.nan
    error_best_distance_minus[use_near_distance] = catalog['error_near_distance_minus'][use_near_distance]
    error_best_distance_minus[use_far_distance] = catalog['error_far_distance_minus'][use_far_distance]
    error_best_distance_minus[ambiguous_distance] = np.nan

    KDA_resolution[use_near_distance] = 'N'
    KDA_resolution[use_far_distance] = 'F'
    KDA_resolution[same_distances] = 'U'
    KDA_resolution[ambiguous_distance] = 'A'

    return best_distance, KDA_resolution, p_near, p_far, error_best_distance_plus, error_best_distance_minus

def assign_distance_columns(catalog, best_distance, KDA_resolution, p_near, p_far, 
                            error_best_distance_plus, error_best_distance_minus):
    """
    Takes the output of the above function and assigns it to columns in the table.

    """

    catalog['distance'] = best_distance
    catalog['error_distance_plus'] = error_best_distance_plus
    catalog['error_distance_minus'] = error_best_distance_minus
    catalog['KDA_resolution'] = KDA_resolution
    catalog['p_near'] = p_near
    catalog['p_far'] = p_far


def assign_distance_columns_trivial(catalog):
    """ Use this in the 2nd and 3rd quadrants. """

    catalog['distance'] = catalog['distance']
    catalog['error_distance_plus'] = catalog['error_distance_plus']
    catalog['error_distance_minus'] = catalog['error_distance_minus']
    catalog['KDA_resolution'] = np.array(['U']*len(catalog['distance']))
    catalog['p_near'] = np.zeros_like(catalog['distance'])
    catalog['p_far'] = np.zeros_like(catalog['distance'])

