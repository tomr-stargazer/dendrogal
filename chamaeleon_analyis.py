"""
An analysis & catalog production of the Chamaeleon region using Boulanger 1998 as a baseline.

"""

from __future__ import division

import numpy as np

import astropy.units as u
import astropy.constants as c

from astrodendro.scatter import Scatter

from demo import chamaeleon_demo, reduce_selection_to_principal_branches
from assign_physical_values import assign_size_mass_alpha_pressure

from hub_selection_manager import reconstruct_selections_from_template

def generate_chamaeleon_catalog():
	d, catalog, x, y = chamaeleon_demo(downsample_factor=1, resample=False, recenter=False, min_npix=10, min_value=0.1, min_delta=0.1)

	catalog['Distance'] = np.ones(len(catalog)) * 160 * u.pc
	catalog['size'], catalog['mass'], catalog['alpha'], catalog['pressure'] = assign_size_mass_alpha_pressure(catalog)
	catalog['filling_factor'] = catalog['area_exact'] / catalog['area_ellipse']

	return d, catalog, x, y

def load_chamaeleon_catalog_from_template():
	filename = 'Chamaeleon_detail'

	d, catalog, header, metadata = chamaeleon_demo(downsample_factor=1, resample=False, recenter=False, min_npix=10, min_value=0.1, min_delta=0.1)

	(template_cube, template_header, new_selection_dictionary, 
		selection_ID_dictionary) = reconstruct_selections_from_template(d, filename)

	catalog['Distance'] = np.ones(len(catalog)) * 160 * u.pc
	catalog['size'], catalog['mass'], catalog['alpha'], catalog['pressure'] = assign_size_mass_alpha_pressure(catalog)
	catalog['filling_factor'] = catalog['area_exact'] / catalog['area_ellipse']

	return d, catalog, header, metadata, template_cube, template_header, new_selection_dictionary, selection_ID_dictionary

