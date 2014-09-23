"""
An analysis & catalog production of the Orion region using Wilson 2005 as a baseline. 

"""

from __future__ import division

import astropy.units as u
import astropy.constants as c

from demo import orion_demo, reduce_selection_to_principal_branches
# def orion_demo(**kwargs):
#     return downsampled_demo('DHT27_Orion_mom.fits', **kwargs)

from assign_physical_values import assign_size_mass_alpha_pressure
from hub_selection_manager import reconstruct_selections_from_template, assign_region_dict_distances


# d, catalog, x, y = orion_demo(downsample_factor=1, resample=False, recenter=False, min_npix=10, min_value=0.1, min_delta=0.1)

filename = 'Orion_detail'

def load_orion_catalog_from_template():

    d, catalog, header, metadata = orion_demo(downsample_factor=1, resample=False, recenter=False, min_npix=10, min_value=0.1, min_delta=0.1)

    (template_cube, template_header, new_selection_dictionary, 
        selection_ID_dictionary) = reconstruct_selections_from_template(d, filename)

    # each sub-region gets its own distance!

    # Southern Filament: 460 pc
    # Mon R2: 832 pc*
    # NGC 2149: 425
    # Northern Filament: 393
    # Orion A: 483*
    # Scissors: 150
    # Orion East: 120
    # Orion B: 420?

    distance_dict = {}
    distance_dict['Southern Filament'] = 460
    distance_dict['Mon R2'] = 832
    distance_dict['Orion A'] = 483
    distance_dict['Orion B'] = 420
    distance_dict['Orion East'] = 120
    distance_dict['NGC 2149'] = 425
    distance_dict['Scissors'] = 150
    distance_dict['Northern Filament'] = 393

    distances = assign_region_dict_distances(catalog, new_selection_dictionary, distance_dict)

    catalog['Distance'] = distances

    catalog['size'], catalog['mass'], catalog['alpha'], catalog['pressure'] = assign_size_mass_alpha_pressure(catalog)
    catalog['filling_factor'] = catalog['area_exact'] / catalog['area_ellipse']

    return d, catalog, header, metadata, template_cube, template_header, new_selection_dictionary, selection_ID_dictionary

