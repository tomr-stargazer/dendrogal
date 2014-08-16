"""
Little functions that help manage dendrogram selections.

"""

from __future__ import division

import numpy as np

import astrodendro

# some of these will wanna be function generators

def transfer_function_generator(hub, origin=1, destination=2):

    def transfer_between_selections():

        hub.select(destination, hub.selections[origin]+hub.selections[destination], subtree=False)
        hub.select(origin, None)

    return transfer_between_selections

selection_dictionary = {}
selection_ID_dictionary = {}

def stash_function_generator(hub, origin=2, destination=3):

    def stash_selection(name):

        if type(name) is not str:
            raise ValueError("`name` must be a string, e.g. 'Orion'")

        try:
            selection_dictionary[name] += hub.selections[origin]
        except KeyError:
            selection_dictionary[name] = hub.selections[origin]
            selection_ID_dictionary[name] = len(selection_dictionary.keys()) # increment by one

        hub.select(destination, hub.selections[origin]+hub.selections[destination], subtree=False)
        hub.select(origin, None)

    return stash_selection

