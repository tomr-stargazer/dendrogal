"""
Quick and dirty reproduction of Rosolowsky et al. 2008's analysis of GMCs in orion/Monoceros

"""

from __future__ import division

import numpy as np
import astropy
import astropy.units as u

from demo import orion_demo
from assign_physical_values import assign_size_mass_alpha_pressure

d, catalog, datacube_dt_header, metadata = orion_demo(downsample_factor=1, min_npix=20)

# Orion structures

# Orion A trunk

oriona_trunk = [structure for structure in d.all_structures if structure.idx == 47][0]
oriona_structures = oriona_trunk.descendants
oriona_indices = [structure.idx for structure in oriona_structures]

# probably the wrong approach
oriona_catalog = catalog[np.in1d(catalog['_idx'], oriona_indices)]

# Orion B trunk

orionb_trunk = [structure for structure in d.all_structures if structure.idx == 58][0]
orionb_structures = orionb_trunk.descendants
orionb_indices = [structure.idx for structure in orionb_structures]

# see above - wrong approach
orionb_catalog = catalog[np.in1d(catalog['_idx'], orionb_indices)]

# Monoceros trunk

monoceros_trunk = [structure for structure in d.all_structures if structure.idx == 158][0]
monoceros_structures = monoceros_trunk.descendants
monoceros_indices = [structure.idx for structure in monoceros_structures]

monoceros_catalog = catalog[np.in1d(catalog['_idx'], monoceros_indices)]

# okay. Let's apply distances and then do virial & pressure analysis

distance_array = np.ones_like(catalog['_idx'])
distance_array[np.in1d(catalog['_idx'], oriona_indices)] = 450
distance_array[np.in1d(catalog['_idx'], orionb_indices)] = 450
distance_array[np.in1d(catalog['_idx'], monoceros_indices)] = 800

distance_column = astropy.table.Column(
	data=distance_array*u.pc,
	name="Distance")
catalog.add_column(distance_column)

s, m, v, p = assign_size_mass_alpha_pressure(catalog)

catalog['size'] = astropy.table.Column(data=s, name='size')
catalog['mass'] = astropy.table.Column(data=m, name='mass')
catalog['virial'] = astropy.table.Column(data=v, name='virial')
catalog['pressure'] = astropy.table.Column(data=p, name='pressure')

#
