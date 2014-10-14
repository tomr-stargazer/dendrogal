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

# the plot

xs = np.logspace(-1,2,20)
ys = myoutput.beta[0] * xs**myoutput.beta[1]

dv_orion = d_orion.viewer()

dv_orion.hub.select(selection)

dsl_orion = Scatter(d_orion, dv_orion.hub, catalog_orion, 'size', 'v_rms')
dsl_orion.set_loglog()

plt.plot(xs, ys, scalex=False, scaley=False)