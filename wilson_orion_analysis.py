"""
An analysis & catalog production of the Orion region using Wilson 2005 as a baseline. 

"""

from __future__ import division

import astropy.units as u
import astropy.constants as c

from demo import orion_demo
# def orion_demo(**kwargs):
#     return downsampled_demo('DHT27_Orion_mom.fits', **kwargs)

d, catalog, x, y = orion_demo(downsample_factor=2, resample=True, recenter=False, min_npix=10, min_value=0.1, min_delta=0.1)

