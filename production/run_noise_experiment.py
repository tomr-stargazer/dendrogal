"""
This is a module that aims to address how our results vary when random noise is added.

It focuses on the first quadrant, as a representative subset of the data.

Basic ideas:

1. run the analysis un-molested and measure the output. 
(Size-linewidth results, mass spectrum results, total mass, number of clouds.)

2. Create a lightly salted version of the data and repeat the above, just once though.

3. Once that's all established, do it like 20-100 times (maybe just ten?) 
   with 0.5x, 1.0x, and 2.0x the noise.


"""

from __future__ import division

import os
import numpy as np

from dendrogal.production.convenience_function import load_permute_dendro_catalog, load_permute_data
from dendrogal.production.load_and_process_data import load_data, permute_data_to_standard_order

from dendrogal.config import data_path

# btw the noise here is 0.18 K / px / channel
data_filename = "DHT08_Quad1_interp.fits"

# so for the "part 2" thing above... we'll have to take the INTERPOLATED data and moment-mask it manually yeah?

def do_a_thing(noise_added=0):

    if noise_added > 0:

        # get some data
        datacube, header = load_permute_data(data_filename)

        # add some noise to it
        noise_cube = np.random.normal(scale=noise_added, size=datacube.shape)
        noise_added_cube = datacube + noise_cube

        # moment-mask it






        # add the 
    else:

        from dendrogal.production.make_firstquad_stub import d, catalog