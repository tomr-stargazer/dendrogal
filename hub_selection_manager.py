"""
Little functions that help manage dendrogram selections.

"""

from __future__ import division
import os

import random

import numpy as np

import astrodendro
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

def stash_function_generator(hub, origin=2, destination=3):

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

    return stash_selection, selection_dictionary, selection_ID_dictionary

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

def save_template_cube(cube, header, clobber=True, filename='template', output_path=None):

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

    print "Template file saved to {0}".format(output_path)

    return output_path
