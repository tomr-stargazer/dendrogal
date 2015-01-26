"""
Visually displays plots and fits from catalog_measurement.py

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from .catalog_measurement import size_linewidth_slope


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


