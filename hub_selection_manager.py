"""
Little functions that help manage dendrogram selections.

"""

from __future__ import division

import random

import numpy as np

import astrodendro

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


