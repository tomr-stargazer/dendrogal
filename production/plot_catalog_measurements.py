"""
Visually displays plots and fits from catalog_measurement.py

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

from dendrogal.production.catalog_measurement import (
    size_linewidth_slope, cumulative_massfunction_fit, 
    truncated_cloudmass_function, powerlaw_cloudmass_function)


def plot_size_linewidth_fit(catalog, title=""):

    size_linewidth_output = size_linewidth_slope(catalog)
    size_array = catalog['size']
    linewidth_array = catalog['v_rms']

    fit_coefficient = size_linewidth_output.beta[0]
    fit_exponent = size_linewidth_output.beta[1]

    fit_xs = np.logspace(-2, 8, 20)
    fit_ys = fit_coefficient * fit_xs ** fit_exponent

    fig = plt.figure()

    plt.plot(size_array, linewidth_array, 'ro')
    plt.loglog()
    plt.plot(fit_xs, fit_ys, 'g--', zorder=0, scalex=False, scaley=False,
             label=r"$\sigma_v$ = {0:.2f} $\times$ R$^{{{1:.2f}}}$".format(fit_coefficient, fit_exponent))

    plt.legend(loc='upper left')

    plt.xlabel("cloud radius R (pc)", fontsize=18)
    plt.ylabel("cloud linewdth $\sigma_v$ (km s$^{-1}$)", fontsize=18)
    if title is not "":
    	plt.title(title)

    return fig, size_linewidth_output

def plot_cmf_without_fit(catalog):

	fig = plt.figure()

	number_in_bin_q2, bin_edges, ch = plt.hist(np.log10(catalog['mass']), cumulative=-1, log=True, bins=20)
	bin_centers_q2 = (bin_edges[1:] + bin_edges[:-1])/2.
	plt.clf()
	plt.plot(bin_centers_q2, number_in_bin_q2, 'ko' )
	plt.semilogy()
	plt.xlabel(r"log$_{10}$ (M$_{GMC}$ / M$_\odot$)")
	plt.ylabel("n(M > M')")

	return fig

def plot_cmf(catalog, bins=20, hist_range=(4, 7), mass_column_name='mass', **kwargs):

    fig = plt.figure()

    number_in_bin_q2, bin_edges, ch = plt.hist(np.log10(catalog[mass_column_name]), cumulative=-1, log=True, bins=bins, range=hist_range)
    bin_centers_q2 = (bin_edges[1:] + bin_edges[:-1])/2.
    plt.clf()
    plt.plot(bin_centers_q2, number_in_bin_q2, 'ko' )
    plt.semilogy()
    plt.xlabel(r"log$_{10}$ (M$_{GMC}$ / M$_\odot$)")
    plt.ylabel("n(M > M')")
    cmf_output = cumulative_massfunction_fit(catalog, bins=bins, mass_column_name=mass_column_name, **kwargs)

    M_0, N_0, gamma = cmf_output[0]

    m_array = np.linspace(min(bin_edges), max(bin_edges), 50)
    n_array = truncated_cloudmass_function([M_0, N_0, gamma], 10**m_array)

    plt.plot(m_array, n_array, label="$\\gamma = {0:.2f}$,\n$M_0={1:.2e}$,\n$N_0={2:.1f}$".format(gamma, M_0, N_0))

    text_string = r"$N(M' > M) = N_0 \left [ \left ( \frac{M}{M_0} \right )^{\gamma+1} - 1 \right ]$"

    plt.text(4.1, 3, text_string, fontsize=18)
    plt.xlim(*hist_range)
    plt.ylim(0.7, 1e3)

    plt.legend(loc='upper right')

    return fig, cmf_output


def plot_cmf_with_arbitrary_input(catalog, cmf_input, bins=20, hist_range=(4, 7), 
                                  mass_column_name='mass', **kwargs):

    fig = plt.figure()

    number_in_bin_q2, bin_edges, ch = plt.hist(np.log10(catalog[mass_column_name]), cumulative=-1, log=True, bins=bins, range=hist_range)
    bin_centers_q2 = (bin_edges[1:] + bin_edges[:-1])/2.
    plt.clf()
    plt.plot(bin_centers_q2, number_in_bin_q2, 'ko' )
    plt.semilogy()
    plt.xlabel(r"log$_{10}$ (M$_{GMC}$ / M$_\odot$)")
    plt.ylabel("n(M > M')")

    M_0, N_0, gamma = cmf_input

    m_array = np.linspace(min(bin_edges), max(bin_edges), 50)
    n_array = truncated_cloudmass_function([M_0, N_0, gamma], 10**m_array)

    plt.plot(m_array, n_array, label="$\\gamma = {0:.2f}$,\n$M_0={1:.2e}$,\n$N_0={2:.1f}$".format(gamma, M_0, N_0))

    text_string = r"$N(M' > M) = N_0 \left [ \left ( \frac{M}{M_0} \right )^{\gamma+1} - 1 \right ]$"

    plt.text(4.1, 3, text_string, fontsize=18)
    plt.xlim(*hist_range)
    plt.ylim(0.7, 1e3)

    plt.legend(loc='upper right')

    return fig    


def plot_mass_density_radius(catalog, bins=20):

    fig = plt.figure()

    catalog_rgal = catalog['R_gal']

    radii = np.arange(bins)
    mass_per_radius = np.zeros_like(radii) 
    area_per_radius = np.zeros_like(radii) * u.kpc**2

    for i, r in enumerate(radii):

        if i == bins-1:
            break

        catalog_at_r = catalog[(catalog_rgal > radii[i]) & 
                               (catalog_rgal < radii[i+1])]

        mass_in_r = np.sum(catalog_at_r['mass'])

        mass_per_radius[i] = mass_in_r

        area_per_radius[i] = np.pi*((radii[i+1])**2 - radii[i]**2) * u.kpc**2

    mass_density_per_radius = (mass_per_radius*u.solMass/area_per_radius).to(u.solMass / u.pc**2)

    plt.plot(radii, mass_density_per_radius, 'ko-')
    plt.ylabel('Mass in each kpc annulus, divided by total annulus area')
    plt.xlabel('Galactocentric distance (kpc)')

    plt.xlim(2,15)

    fig2 = plt.figure()

    plt.plot(radii, mass_per_radius, 'ko-')
    plt.semilogy()
    plt.ylabel('Total mass per kpc annulus')
    plt.xlabel('Galactocentric distance (kpc)')

    plt.xlim(2,15)

    return mass_per_radius, mass_density_per_radius


def plot_size_linewidth_with_nearfar(catalog, ax, labels=True, distance_threshold=8, alternate_style=False):

    unambiguous_nearby = (catalog['KDA_resolution'] == 'U') & (catalog['distance'] < distance_threshold)
    unambiguous_far = (catalog['KDA_resolution'] == 'U') & (catalog['distance'] >= distance_threshold)
    near = catalog['KDA_resolution'] == 'N'
    far = catalog['KDA_resolution'] == 'F'

    size_linewidth_output = size_linewidth_slope(catalog[near|far|unambiguous_nearby])

    size_array = catalog['size']
    linewidth_array = catalog['v_rms']

    fit_coefficient = size_linewidth_output.beta[0]
    sd_coefficient = size_linewidth_output.sd_beta[0]
    fit_exponent = size_linewidth_output.beta[1]
    sd_exponent = size_linewidth_output.sd_beta[1]

    fit_xs = np.logspace(-2, 8, 20)
    fit_ys = fit_coefficient * fit_xs ** fit_exponent

    if alternate_style:
        label_string = 'tangent'
    else:
        label_string = 'unambiguous (nearby)'

    ax.plot(np.log10(size_array[unambiguous_nearby]), 
            np.log10(linewidth_array[unambiguous_nearby]), 'k.', 
            ms=6, label=label_string)

    if alternate_style:
        pass
    else:
        ax.plot(np.log10(size_array[unambiguous_far]), np.log10(linewidth_array[unambiguous_far]), 'wo', label='unambiguous (far)', ms=3)

    ax.plot(np.log10(size_array[near]), np.log10(linewidth_array[near]), 'o', markerfacecolor='none', markeredgecolor='b', mew=0.5, ms=3.5, label='near')
    ax.plot(np.log10(size_array[far]), np.log10(linewidth_array[far]), 'o', markerfacecolor='none', markeredgecolor='r', mew=0.5, ms=3.5, label='far')

    if labels:
        fit_label = r"$\sigma_v$ = {0:.2f}$\pm${2:.2f} $\times$ R$^{{{1:.2f} \pm {3:.2f}}}$".format(fit_coefficient, fit_exponent, sd_coefficient, sd_exponent) 
    elif alternate_style:
        fit_label = '_nolegend_'
    else:
        fit_label = r"$\sigma_v = A \times R^\beta$"

    ax.plot(np.log10(fit_xs), np.log10(fit_ys), 'g--', zorder=0, scalex=False, scaley=False,
             label=fit_label)

    if labels:
        ax.legend(loc='lower right')

        ax.set_xlabel(r"$\log (R/\rm{pc})$", fontsize=18)
        ax.set_ylabel(r"$\log (\sigma_v / \rm{km s}^{-1})$", fontsize=18)

    return size_linewidth_output

def plot_size_linewidth_with_nearfar_fig(catalog):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    return fig, plot_size_linewidth_with_nearfar(catalog, ax)


def plot_cmf_except_farU(catalog, ax, bins=20, hist_range=(4, 7), labels=True, **kwargs):

    unambiguous_far = (catalog['KDA_resolution'] == 'U') & (catalog['distance'] >= 8)

    number_in_bin_q2, bin_edges, ch = ax.hist(np.log10(catalog[~unambiguous_far]['mass']), cumulative=-1, log=True, bins=bins, range=hist_range)
    bin_centers_q2 = (bin_edges[1:] + bin_edges[:-1])/2.
    # plt.clf()
    ax.plot(bin_centers_q2, number_in_bin_q2, 'ko' )
    ax.semilogy()
    if labels:
        ax.set_xlabel(r"log$_{10}$ (M$_{GMC}$ / M$_\odot$)")
        ax.set_ylabel("n(M > M')")

    cmf_output = cumulative_massfunction_fit(catalog[~unambiguous_far], bins=bins, mass_column_name='mass', **kwargs)

    M_0, N_0, gamma = cmf_output[0]

    m_array = np.linspace(min(bin_edges), max(bin_edges), 50)
    n_array = truncated_cloudmass_function([M_0, N_0, gamma], 10**m_array)

    plt.plot(m_array, n_array, label="$\\gamma = {0:.2f}$,\n$M_0={1:.2e}$,\n$N_0={2:.1f}$".format(gamma, M_0, N_0))

    if labels:
        text_string = r"$N(M' > M) = N_0 \left [ \left ( \frac{M}{M_0} \right )^{\gamma+1} - 1 \right ]$"

        ax.text(4.1, 3, text_string, fontsize=18)
        ax.legend(loc='upper right')

    ax.set_xlim(*hist_range)
    ax.set_ylim(0.7, 1e3)


    return cmf_output


def plot_cmf_with_fit_fig(catalog, **kwargs):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    return fig, plot_cmf_except_farU(catalog, ax, **kwargs)


def plot_cmf_with_pl_and_tpl(catalog, cmf_input, ax, bins=20, hist_range=(4, 7), labels=True, **kwargs):

    N_0, tpl_M_0, tpl_gamma, pl_M_0, pl_gamma = cmf_input

    number_in_bin_q2, bin_edges, ch = ax.hist(np.log10(catalog['mass']), cumulative=-1, log=True, bins=bins, range=hist_range)
    bin_centers_q2 = (bin_edges[1:] + bin_edges[:-1])/2.

    ax.plot(bin_centers_q2, number_in_bin_q2, 'ko', ms=3 )
    ax.semilogy()
    # if labels:
    #     ax.set_xlabel(r"log$_{10}$ (M$_{GMC}$ / M$_\odot$)")
    #     ax.set_ylabel("n(M > M')")

    m_array = np.linspace(min(bin_edges), max(bin_edges), 50)
    tpl_n_array = truncated_cloudmass_function([tpl_M_0, N_0, tpl_gamma], 10**m_array)
    pl_n_array = powerlaw_cloudmass_function([pl_M_0, pl_gamma], 10**m_array)

    ax.plot(m_array, tpl_n_array, 'b-', lw=2, label="TPL $\\gamma = {0:.2f}$,\n$M_0={1:.2e}$,\n$N_0={2:.1f}$".format(tpl_gamma, tpl_M_0, N_0))
    ax.plot(m_array, pl_n_array, 'r--', lw=2, label="PL $\\gamma = {0:.2f}$,\n$M_0={1:.2e}$".format(pl_gamma, pl_M_0))

    if labels:
        ax.legend(loc='upper right')

    ax.set_xlim(*hist_range)
    ax.set_ylim(0.7, 1e3)

    return 


def plot_cmf_with_pl_and_tpl_fig(catalog, cmf_input, **kwargs):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    return fig, plot_cmf_with_pl_and_tpl(catalog, cmf_input, ax, **kwargs)

