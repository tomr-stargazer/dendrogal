"""
This is a special version of `demo.py` that contains 
only things Matt P needs for his Astro 98 project.

"""

# first, let's make sure that everything imports properly

from __future__ import division

import os.path
import sys
from functools import partial
import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs
from astropy.io.fits import getdata, getheader
import astropy.io.fits as fits

matt_path = os.path.expanduser("~/Desktop/PythonWorkspace/dendrogal/")
matt_code_path = os.path.expanduser("~/Desktop/PythonWorkspace/")

try:
    from dendrogal.demo import downsampled_demo
except ImportError:
    sys.path.append(matt_code_path)
    from dendrogal.demo import downsampled_demo

from dendrogal.integrated_and_cartoon_viewer import CartoonDualViewer
import dendrogal.hub_selection_manager as selection_manager
from dendrogal.hub_selection_manager import (
  transfer_function_generator, stash_function_generator, assign_region_dict_distances,
  save_template, load_template)
from dendrogal.draw_schlafly_points_on_map import draw_big_cloud_table, draw_mbm_cloud_table


matt_demo_function = partial(
    downsampled_demo, 
    data_file='COGAL_local_mom.fits', 
    downsample_factor=2,
    min_npix=20,
    min_delta=0.1,
    min_value=0.1,
    data_path=matt_path,
    memmap=False)

def run_matt_demo():

    d, catalog, header, metadata = matt_demo_function()

    dv = d.viewer()

    cdv = CartoonDualViewer(d, dv.hub)

    t = transfer_function_generator(dv.hub)
    s, selection_dictionary, selection_ID_dictionary = stash_function_generator(dv.hub, checkpoint=False)

    output_dict = {}
    output_dict['dendrogram'] = d
    output_dict['catalog'] = catalog
    output_dict['header'] = header
    output_dict['viewer'] = dv
    output_dict['dual_viewer'] = cdv
    output_dict['transfer_function'] = t
    output_dict['stash_function'] = s
    output_dict['selection_dictionary'] = selection_dictionary
    output_dict['selection_ID_dictionary'] = selection_ID_dictionary

    print next_string

    return output_dict


def save_matt_template(output_dict, filename):
    if type(filename) is not str:
        raise ValueError("`filename` argument must be a string")
    if type(output_dict) is not dict:
        raise ValueError("`output_dict` argument must be a dict")

    output_path = matt_path+filename

    save_template(output_dict['selection_dictionary'], output_dict['selection_ID_dictionary'], output_dict['header'], output_path=output_path)


def load_matt_template(output_dict, filename, dv=None):

    input_path = matt_path+filename

    (cube, header, selection_idx_dictionary, 
        selection_ID_dictionary) = load_template(filename, input_path, dv=dv)

    output_dict['selection_dictionary'] = selection_idx_dictionary
    output_dict['selection_ID_dictionary'] = selection_ID_dictionary


print(
"""


  Hi Matt - here are some ideas for how to use this tool.
  Start by running the following line:

.. output = run_matt_demo()

  If that works, we are in good shape. Let it run for a minute or so. 
  It will print out some status updates while it computes a dendrogram.""")

next_string = """


  It should load up a couple of interactive graphical panels.
  Let's use the one with the parallel "cartoon" and "data" views 
  of the Galaxy.

  Next, what you should do is initialize your 
  "t" (transfer) and
  "s" (stash or save) 
  functions like this:

.. t = output['transfer_function']
.. s = output['stash_function']

  From this point on, you can select regions in RED in the cartoon-viewer
  (by left-clicking on them),
  then transfer them to your BLUE selection by running the 
  following command

.. t()

  and finally, you can "stash" your BLUE selection into a local 
  "dictionary" data structure (thus making it GREEN) by saying

.. s("Cloud Name Goes Here") # for example, you could say "Orion A"

  You can then access that "dictionary" data structure as

.. selection_dictionary = output['selection_dictionary']

  Let's give this a shot and see where it goes!"""


def draw_big_clouds(output):
    viewer = output['dual_viewer']

    return draw_big_cloud_table(viewer.ax_cartoon)


def draw_mbm_clouds(output):
    viewer = output['dual_viewer']

    return draw_mbm_cloud_table(viewer.ax_cartoon)



