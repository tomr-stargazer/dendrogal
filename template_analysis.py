"""
template_analysis.py : a module for analyzing stuff based on templates.

"""

from __future__ import division

from demo import cogal_local_resampled_demo, cogal_deep_resampled_demo, reduce_selection_to_principal_branches
from hub_selection_manager import reconstruct_selections_from_template

def template_analysis_local_test():

	# let's do the local cube
	filename = 'checkpoint/22:Gem OB1'

    # intentionally SLIGHTLY altered from the original parameters
	# d, catalog, header, metadata = cogal_local_resampled_demo(2, min_npix=10, min_value=0.1, min_delta=0.1)
	d, catalog, header, metadata = cogal_local_resampled_demo(2, min_npix=15, min_value=0.09, min_delta=0.11)

	(template_cube, template_header, new_selection_dictionary, 
		selection_ID_dictionary) = reconstruct_selections_from_template(d, filename)

	return d, catalog, header, metadata, template_cube, template_header, new_selection_dictionary, selection_ID_dictionary

def template_analysis_deep_test():

	filename = 'checkpoint/2:Outer Arm (North)'

    # intentionally SLIGHTLY altered from the original parameters
	# d, catalog, header, metadata = cogal_local_resampled_demo(2, min_npix=10, min_value=0.1, min_delta=0.1)
	d, catalog, header, metadata = cogal_deep_resampled_demo(2, min_npix=15, min_value=0.09, min_delta=0.11)

	(template_cube, template_header, new_selection_dictionary, 
		selection_ID_dictionary) = reconstruct_selections_from_template(d, filename)

	return d, catalog, header, metadata, template_cube, template_header, new_selection_dictionary, selection_ID_dictionary


def traverse_recursively(catalog, structure, depth, function):

	tabs = '.'*depth
	me = "{1}{0}\n".format(function(structure, catalog), tabs)

	if len(structure.children) == 2:
		try:
			return me + traverse_recursively(catalog, structure.children[0], depth+1, function) + traverse_recursively(catalog, structure.children[1], depth+1, function)
		except IndexError as e:
			print "  glitch at {0}: {1}  ".format(structure.idx, e)
			return me
	else:
		return me


def struct_printer(struct, catalog):

	row = catalog['_idx'] == struct.idx

	l = catalog['x_cen'][row][0]
	b = catalog['y_cen'][row][0]
	v = catalog['v_cen'][row][0]

	size = catalog['radius'][row][0]
	distance = catalog['Distance'][row][0]
	mass = catalog['mass'][row][0]

	return "name: {0}, l: {1:.2f}, b: {2:.2f}, v: {3:.2f} km/s,  size: {4:.2f} deg, distance: {5:.2f} pc, mass: {6:.1e} Msun".format(struct.idx, l, b, v, size, distance, mass)

def traverse_and_print_catalog(catalog, columns, structures):

	principal_branches = reduce_selection_to_principal_branches(structures)

	function = struct_printer

	for struct in principal_branches:

		print traverse_recursively(catalog, struct, 0, struct_printer)



