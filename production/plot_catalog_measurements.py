"""
Visually displays plots and fits from catalog_measurement.py

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

from astrodendro_analysis.production.catalog_measurement import size_linewidth_slope, cumulative_massfunction_fit, truncated_cloudmass_function


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

def plot_mass_density_radius(catalog, bins=20):

    fig = plt.figure()

    catalog_rgal = (catalog['x_gal']**2 + catalog['y_gal']**2)**(1/2)

    radii = np.arange(bins)
    mass_per_radius = np.zeros_like(radii) 
    area_per_radius = np.zeros_like(radii)

    for i, r in enumerate(radii):

        if i == bins-1:
            break

        catalog_at_r = catalog[(catalog_rgal > radii[i]) & 
                               (catalog_rgal < radii[i+1])]

        mass_in_r = np.sum(catalog_at_r['mass'])

        mass_per_radius[i] = mass_in_r

        area_per_radius = np.pi*((radii[i+1])**2 - radii[i]**2) * u.kpc**2

    mass_density_per_radius = (mass_per_radius*u.solMass/area_per_radius).to(u.solMass / u.pc**2)

    plt.plot(radii, mass_density_per_radius, 'ko-')

    plt.xlim(2,17)

    return mass_per_radius, mass_density_per_radius
