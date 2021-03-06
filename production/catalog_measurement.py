"""
Functions to measure overall cloud statistics, given a sanitized input cloud catalog.

"""

from __future__ import division

import numpy as np

from scipy.odr import RealData, Model, ODR
from scipy.optimize import leastsq


def powerlaw_function(B, x):
    return B[0] * (x**B[1])


def size_linewidth_slope(catalog):
    """
    Measures a catalog's slope of sizes and linewidths using scipy.odr.

    ODR means "orthogonal distance regression", and is a way of fitting
    models to data where both the x and y values have scatter.

    Parameters
    ----------
    catalog : astropy.table.table.Table
        Table containing columns for 'size' and 'v_rms'.

    Returns
    -------
    odr_output : scipy.odr.odrpack.Output
        Output of the scipy ODR routine. The fit parameters
        are stored in odr_output.beta .

    """

    if 'size' not in catalog.colnames or 'v_rms' not in catalog.colnames:
        raise ValueError("'size' and 'v_rms' must be columns in `catalog`!")

    if 'error_size_plus' not in catalog.colnames or 'error_size_minus' not in catalog.colnames:
        mean_size_error = None
    else:
        mean_size_error = (catalog['error_size_plus']*catalog['error_size_minus'])**(1/2)
        # if any clouds have no error in their size, the procedure will puke unless we do this:
        mean_size_error[mean_size_error==0] = 0.01

    size_linewidth_data = RealData(catalog['size'], catalog['v_rms'], sx=mean_size_error)

    powerlaw_model = Model(powerlaw_function)
    # The provided `beta0` is a throwaway guess that shouldn't change the outcome.
    odr_object = ODR(size_linewidth_data, powerlaw_model, beta0=[1., 1.])
    odr_output = odr_object.run()

    return odr_output


def truncated_cloudmass_function(parameter_list, mass_array):
    """
    The cumulative mass function from Columbo et al. (2014),
    equation 16.

    Parameters
    ----------
    parameter_list : list of floats
        contains M_0, N_0, gamma
        M_0 : the maximumum mass in the distribution
        N_0 : the number of clouds more massive than 2**(1/(gamma+1))*M_0
        gamma : the index describing how mass is distributed

    mass_array : array of x values
        represents the mass bins

    Returns
    -------
    N(M' > M) : array
        number of clouds expected in each mass bin

    """

    M_0, N_0, gamma = parameter_list

    N_by_mass = N_0 * ((mass_array/M_0)**(gamma+1) - 1)

    return N_by_mass


def powerlaw_cloudmass_function(parameter_list, mass_array):
    """
    The cumulative mass function from Columbo et al. (2014),
    equation 15.

    Parameters
    ----------
    parameter_list : list of floats
        contains M_0, gamma
        M_0 : the maximumum mass in the distribution
        gamma : the index describing how mass is distributed

    mass_array : array of x values
        represents the mass bins

    Returns
    -------
    N(M' > M) : array
        number of clouds expected in each mass bin

    """

    M_0, gamma = parameter_list

    N_by_mass = ((mass_array/M_0)**(gamma+1) )

    return N_by_mass


def cumulative_massfunction_fit(catalog, min_mass=1e5, max_mass=3e7, bins=20, mass_column_name='mass'):
    """
    Measures a catalog's mass function slope using .

    ODR means "orthogonal distance regression", and is a way of fitting
    models to data where both the x and y values have scatter.

    Parameters
    ----------
    catalog : astropy.table.table.Table
        Table containing columns for 'size' and 'v_rms'.

    Returns
    -------
    odr_output : scipy.odr.odrpack.Output
        Output of the scipy ODR routine. The fit parameters
        are stored in odr_output.beta .

    """

    if 'mass' not in catalog.colnames:
        raise ValueError("'mass' must be column in `catalog`!")

    # make a reversely-summed histogram (the [::-1]'s ensure it's reversed properly)
    hist, bin_edges = np.histogram( np.log10(catalog[mass_column_name]), range=[np.log10(min_mass), np.log10(max_mass)], bins=bins)
    cum_hist = np.cumsum(hist[::-1])[::-1]
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2.

    initial_gamma = -2
    initial_M_0 = sorted(catalog[mass_column_name])[-2] # skip the biggest, go for second biggest
    initial_N_0 = len(np.where(catalog[mass_column_name] > 2**(1/(initial_gamma+1))*initial_M_0)[0])

    print "initial guess gamma: {0}".format(initial_gamma)
    print "initial guess M_0: {0}".format(initial_M_0)
    print "initial guess N_0: {0}".format(initial_N_0)

    initial_parameters = [initial_M_0, initial_N_0, initial_gamma]

    errfunc = lambda p, x, y: truncated_cloudmass_function(p, x) - y # Distance to the target function

    fit = leastsq(errfunc, initial_parameters, args=(10**bin_centers, cum_hist), full_output=1)

    return fit








