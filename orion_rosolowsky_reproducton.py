"""
Quick and dirty reproduction of Rosolowsky et al. 2008's analysis of GMCs in orion/Monoceros

"""

from __future__ import division

import numpy as np
import astropy
import astropy.units as u

from demo import orion_demo
from assign_physical_values import assign_size_mass_alpha

d, catalog, datacube_dt_header, metadata = orion_demo(downsample_factor=1, min_npix=20)

# Orion structures

# Orion A trunk

oriona_trunk = [structure for structure in d.all_structures if structure.idx == 47][0]
oriona_structures = oriona_trunk.descendants
oriona_indices = [structure.idx for structure in oriona_structures]

oriona_catalog = catalog[np.in1d(catalog['_idx'], oriona_indices)]

# Orion B trunk

orionb_trunk = [structure for structure in d.all_structures if structure.idx == 58][0]
orionb_structures = orionb_trunk.descendants
orionb_indices = [structure.idx for structure in orionb_structures]

orionb_catalog = catalog[np.in1d(catalog['_idx'], orionb_indices)]

# Monoceros trunk

monoceros_trunk = [structure for structure in d.all_structures if structure.idx == 158][0]
monoceros_structures = monoceros_trunk.descendants
monoceros_indices = [structure.idx for structure in monoceros_structures]

monoceros_catalog = catalog[np.in1d(catalog['_idx'], monoceros_indices)]

# okay. Let's apply distances and then do virial & pressure analysis

oriona_distance_column = astropy.table.Column(
	data=np.ones_like(oriona_indices)*u.pc*450, 
	name="Distance")
oriona_catalog.add_column(oriona_distance_column)



orionb_distance_column = astropy.table.Column(
	data=np.ones_like(orionb_indices)*u.pc*450, 
	name="Distance")
orionb_catalog.add_column(orionb_distance_column)

monoceros_distance_column = astropy.table.Column(
	data=np.ones_like(monoceros_indices)*u.pc*800, 
	name="Distance")
monoceros_catalog.add_column(monoceros_distance_column)


for catalog in [oriona_catalog, orionb_catalog, monoceros_catalog]:

	s, m, v = assign_size_mass_alpha(catalog)

	catalog['size'] = astropy.table.Column(data=s, name='size')
	catalog['mass'] = astropy.table.Column(data=m, name='mass')
	catalog['virial'] = astropy.table.Column(data=v, name='virial')

	catalog['pressure'] = catalog['mass'] * catalog['v_rms']**2 / catalog['size']**3

#
