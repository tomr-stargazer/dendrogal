""" 
Script to generate relevant things for my second-year talk.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import RealData, Model, ODR

from astrodendro.scatter import Scatter

# Part 1. Orion complex: 3 GMCs

from wilson_orion_analysis import load_orion_catalog_from_template
from second_quadrant_analysis import secondquad_distance_demo, thirdquad_distance_demo
from demo import reduce_selection_to_principal_branches

orion_output = load_orion_catalog_from_template()
d_orion, catalog_orion, header, metadata, template_cube, template_header, new_selection_dictionary, selection_ID_dictionary = orion_output

# the fit

selection = new_selection_dictionary['Mon R2'] + new_selection_dictionary['Orion A'] + new_selection_dictionary['Orion B']
selection_ids = [struct.idx for struct in selection]

s_size = catalog_orion[selection_ids]['size']
s_linewidth = catalog_orion[selection_ids]['v_rms']

mydata = RealData(s_size, s_linewidth)
def powerlaw(B, x):
    return B[0] * (x**B[1])
powerlaw_model = Model(powerlaw)
myodr = ODR(mydata, powerlaw_model, beta0=[1.,1.])
myoutput = myodr.run()

orion_fit_constant = myoutput.beta[0]
orion_fit_power = myoutput.beta[1]

# the plot

xs = np.logspace(-1,4,20)
ys = orion_fit_constant * xs**orion_fit_power

dv_orion = d_orion.viewer()

dv_orion.hub.select(1, selection)

dsl_orion = Scatter(d_orion, dv_orion.hub, catalog_orion, 'size', 'v_rms')
dsl_orion.set_loglog()

plt.plot(xs, ys, scalex=False, scaley=False)


# Part 2. Second Quadrant

d2, catalog2, header2, metadata2 = secondquad_distance_demo()

dv_2 = d2.viewer()

selection2_ids = ( (catalog2['v_rms'] < 2*orion_fit_constant * catalog2['size']**orion_fit_power) &
    (catalog2['v_rms'] > 1/3 * orion_fit_constant * catalog2['size']**orion_fit_power) &
    (catalog2['mass'] > 5e4) )
selection2 = [d2[idx] for idx in catalog2['_idx'][selection2_ids]]

gmc_selection2 = reduce_selection_to_principal_branches(selection2)
gmc2_ids = [struct.idx for struct in gmc_selection2]

# plot the cumulative mass function, the size-linewidth dist.

secondquad_data = RealData(catalog2['size'][gmc2_ids], catalog2['v_rms'][gmc2_ids])
secondquad_odr = ODR(secondquad_data, powerlaw_model, beta0=[orion_fit_constant,orion_fit_power]) # use previous results as first guess
secondquad_output = secondquad_odr.run()

quad2_fit_constant = secondquad_output.beta[0]
quad2_fit_power = secondquad_output.beta[1]

xs = np.logspace(-1,4,20)
ys_2 = quad2_fit_constant * xs**quad2_fit_power


plt.figure()
plt.plot(catalog2['size'][gmc2_ids], catalog2['v_rms'][gmc2_ids], 'ro')
plt.plot(xs, ys, '--')

plt.figure()
number_in_bin_q2, bin_edges, ch = plt.hist(np.log10(catalog2['mass'][gmc2_ids]), cumulative=-1, log=True, bins=20)
bin_centers_q2 = (bin_edges[1:] + bin_edges[:-1])/2.
plt.clf()
plt.plot(bin_centers_q2, number_in_bin_q2, 'ko' )
plt.semilogy()
plt.xlabel(r"log$_{10}$ (M$_{GMC}$ / M$_\odot$)")
plt.ylabel("n(M > M')")


# Part 3. Third Quadrant

d3, catalog3, header3, metadata3 = thirdquad_distance_demo()

dv_3 = d3.viewer()

selection3_ids = ( (catalog3['v_rms'] < 2*orion_fit_constant * catalog3['size']**orion_fit_power) &
    (catalog3['v_rms'] > 1/3 * orion_fit_constant * catalog3['size']**orion_fit_power) &
    (catalog3['mass'] > 5e4) )
selection3 = [d3[idx] for idx in catalog3['_idx'][selection3_ids]]

gmc_selection3 = reduce_selection_to_principal_branches(selection3)
gmc3_ids = [struct.idx for struct in gmc_selection3]

# plot the cumulative mass function, the size-linewidth dist.

thirdquad_data = RealData(catalog3['size'][gmc3_ids], catalog3['v_rms'][gmc3_ids])
thirdquad_odr = ODR(thirdquad_data, powerlaw_model, beta0=[orion_fit_constant,orion_fit_power]) # use previous results as first guess
thirdquad_output = thirdquad_odr.run()

quad3_fit_constant = thirdquad_output.beta[0]
quad3_fit_power = thirdquad_output.beta[1]

xs = np.logspace(-1,4,30)
ys_3 = quad3_fit_constant * xs**quad3_fit_power


plt.figure()
plt.plot(catalog3['size'][gmc3_ids], catalog3['v_rms'][gmc3_ids], 'ro')
plt.plot(xs, ys, '--')

plt.figure()
number_in_bin_q3, bin_edges, ch = plt.hist(np.log10(catalog3['mass'][gmc3_ids]), cumulative=-1, log=True, bins=20)
bin_centers_q3 = (bin_edges[1:] + bin_edges[:-1])/2.
plt.clf()
plt.plot(bin_centers_q3, number_in_bin_q3, 'ko' )
plt.semilogy()
plt.xlabel(r"log$_{10}$ (M$_{GMC}$ / M$_\odot$)")
plt.ylabel("n(M > M')")