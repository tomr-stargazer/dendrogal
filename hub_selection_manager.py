"""
Little functions that help manage dendrogram selections.

"""

from __future__ import division
import os
import pickle

import random

import numpy as np
import astropy.units as u

import astrodendro
from astropy.io.fits import getdata, getheader
import astropy.io.fits as fits

dropbox_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/")

# some of these will wanna be function generators

def transfer_function_generator(hub, origin=1, destination=2):

    def transfer_between_selections():

        try:
            joint_selection = set([x for x in hub.selections[origin]+hub.selections[destination] if x is not None])
            hub.select(destination, list(joint_selection), subtree=False)
        except KeyError:
            hub.select(destination, hub.selections[origin], subtree=False)            
        hub.select(origin, None)

    return transfer_between_selections

def stash_function_generator(hub, origin=2, destination=3, checkpoint=True):

    selection_dictionary = {}
    selection_ID_dictionary = {}

    def stash_selection(name):

        if type(name) is not str:
            raise ValueError("`name` must be a string, e.g. 'Orion'")

        try:
            selection_dictionary[name] += hub.selections[origin]
        except KeyError:
            selection_dictionary[name] = hub.selections[origin]
            selection_ID_dictionary[name] = len(selection_dictionary.keys()) # increment by one

        try:
            joint_selection = [x for x in hub.selections[origin]+hub.selections[destination] if x is not None]
            hub.select(destination, joint_selection, subtree=False)
        except KeyError:
            hub.select(destination, hub.selections[origin], subtree=False)            
        hub.select(origin, None)

        # checkpoint system... slow but saves my butt
        if checkpoint:
            filename="{2}/{0}:{1}".format(selection_ID_dictionary[name], name, 'checkpoint')
            save_template_cube(make_template_cube(selection_dictionary, selection_ID_dictionary), None, filename=filename, verbose=False)
            save_template_information(selection_dictionary, selection_ID_dictionary, filename=filename)

    return stash_selection, selection_dictionary, selection_ID_dictionary

def intersection_function_generator(hub, input_a=1, input_b=2, destination=3):

    def intersection_of_selections():

        intersected_selection = set([x for x in hub.selections[input_a] if x is not None]).intersection(set([x for x in hub.selections[input_b] if x is not None]))
        hub.select(destination, list(intersected_selection), subtree=False)

    return intersection_of_selections

def restore_selections_from_dictionary(hub, selection_dictionary, destination=3):
    """ Useful if you mis-click and ruin your green selection. """

    # build a thing from all of the previous selections
    selection_list = []
    for key, value in selection_dictionary.items():
        selection_list.extend(value)

    hub.select(destination, selection_list, subtree=False)

# Now I need "data dump" functions.

def make_template_cube(selection_dictionary, selection_ID_dictionary):

    if len(selection_dictionary.keys()) < 1:
        raise ValueError("Empty selection dictionary!")

    cube = np.zeros_like(selection_dictionary[random.choice(selection_dictionary.keys())][0].get_mask(), dtype=int)

    for key in selection_dictionary.keys():

        print key

        selection = selection_dictionary[key]
        ID = selection_ID_dictionary[key]

        key_cube = np.zeros_like(cube)

        for struct in selection:

            key_cube[struct.indices(subtree=True)] = ID

        cube += key_cube

    return cube

def save_template_cube(cube, header, clobber=True, filename='template', output_path=None, verbose=True):

    output_path = output_path or "{0}templates/{1}.fits".format(dropbox_path, filename)

    try:
        fits.writeto(output_path, cube, header)
    except IOError, e:
        if clobber:
            os.remove(output_path)
            fits.writeto(output_path, cube, header)
        else:
            print "File not saved: {0}".format(e)
            return cube, header

    if verbose:
        print "Template file saved to {0}".format(output_path)

    return output_path

def save_template_information(selection_dictionary, selection_ID_dictionary, filename='pickled_dict', output_path=None):

    output_path = output_path or "{0}templates/{1}".format(dropbox_path, filename)

    selection_idx_dictionary = {}
    for key in selection_dictionary.keys():
        selection_idx_dictionary[key] = [struct.idx for struct in selection_dictionary[key]]

    with open(output_path+"_idx", 'w') as handle:
        pickle.dump(selection_idx_dictionary, handle)

    with open(output_path+"_ID", 'w') as handle:
        pickle.dump(selection_ID_dictionary, handle)

def save_template(selection_dictionary, selection_ID_dictionary, header, **kwargs):
    """ Convenience function - combines `make_template_cube`, `save_template_cube`, and `save_template_information` """

    cube = make_template_cube(selection_dictionary, selection_ID_dictionary)

    save_template_cube(cube, header, **kwargs)
    save_template_information(selection_dictionary, selection_ID_dictionary, **kwargs)

def load_template_cube(filename='template', input_path=None):
    """ Inverts save_template_cube """

    input_path = input_path or "{0}templates/{1}.fits".format(dropbox_path, filename)

    cube, header = getdata(input_path, memmap=False, header=True)

    return cube, header

def _load_template_information(filename='pickled_dict', input_path=None):
    """ Inverts save_template_information """

    input_path = input_path or "{0}templates/{1}".format(dropbox_path, filename)

    with open(input_path+"_idx", 'r') as handle:
        selection_idx_dictionary = pickle.load(handle)

    with open(input_path+"_ID", 'r') as handle:
        selection_ID_dictionary = pickle.load(handle)

    return selection_idx_dictionary, selection_ID_dictionary

def load_template(filename, input_path=None, dv=None, selection_key=3, verbose=True):
    """ Combines _load_template_information and load_template_cube """

    cube, header = load_template_cube(filename=filename, input_path=input_path)
    selection_idx_dictionary, selection_ID_dictionary = _load_template_information(filename=filename, input_path=input_path)

    if dv is not None:
        hub = dv.hub
        dendrogram = dv.dendrogram

        # Beware: no checks are performed to ensure that `dendrogram`
        # is the original dendrogram the template was compiled from

        selection_dictionary = {}
        # build the selection from structures
        for key, idx_list in selection_idx_dictionary.items():

            selection = []
            for idx in idx_list:
                selection.append(dendrogram[idx])
            # selection = dendrogram[idx]
            selection_dictionary[key] = selection

            # select all of these guys
            try:
                joint_selection = set([x for x in selection+hub.selections[selection_key] if x is not None])
                hub.select(selection_key, list(joint_selection), subtree=False)
            except KeyError:
                hub.select(selection_key, selection, subtree=False)            

        return cube, header, selection_dictionary, selection_ID_dictionary

    else:
        if verbose:
            print ("If you provide a `dv` keyword, the template can be "
                   "loaded into the current viewer")

        return cube, header, selection_idx_dictionary, selection_ID_dictionary

# from stack overflow: stackoverflow.com/questions/16216078/test-for-membership-in-a-2d-numpy-array
def inNd(a, b, **kwargs):
    """ Test whether each element of an N-D array is also present in a second array. """
    a = np.asarray(a, order='C')
    b = np.asarray(b, order='C')

    a = a.ravel().view((np.str, a.itemsize * a.shape[1]))
    b = b.ravel().view((np.str, b.itemsize * b.shape[1]))

    return np.in1d(a, b, **kwargs)

def inNd_indices(a, b, **kwargs):
    """ Properly feeds `struct.indices` and `np.where` output into `inNd()` """

    a = np.vstack(a).T
    b = np.vstack(b).T

    return inNd(a, b, **kwargs)

def reconstruct_selections_from_template(dendrogram, filename, input_path=None, **kwargs):
    """ Loads a template, compares it against a dendrogram / dataset, and then assigns structures to template regions """

    # currently assumes input data dimensions match template data dimensions

    (template_cube, template_header, 
        selection_idx_dictionary, 
        selection_ID_dictionary) = load_template(filename=filename, input_path=input_path, verbose=False, **kwargs)

    new_selection_dictionary = {} # name -> structures

    # for each template region...
    for key, value in selection_ID_dictionary.items():

        print "solving for {0}".format(key)
        region_indices = np.where(template_cube == value)

        new_selection_dictionary[key] = []

        for i in range(len(dendrogram)):

            struct = dendrogram[i]
            struct_indices = struct.indices(subtree=True)

            overlap_of_region = inNd_indices(region_indices, struct_indices)

            if overlap_of_region.any():

                # check if MOST of the structure is in the region
                overlap_of_struct = inNd_indices(struct_indices, region_indices)

                if np.sum(overlap_of_struct) > 0.9*len(overlap_of_struct):
                    # add it to the dict
                    new_selection_dictionary[key].append(struct)

    return template_cube, template_header, new_selection_dictionary, selection_ID_dictionary


def assign_region_dict_distances(catalog, selection_dictionary, distance_dictionary):
    """ Assigns distances on a region-by-region basis in a catalog. """

    # assuming that the two dicts have comparable keys:

    distances = np.ones(len(catalog)) * u.pc * np.nan

    for key in selection_dictionary.keys():

        # turn the structures into IDXs and mark them in the catalog
        structures = selection_dictionary[key]
        distance = distance_dictionary[key]
        IDXs = np.array([structure.idx for structure in structures])

        distances[IDXs] = u.Quantity(distance, u.pc)

    return distances


