"""
This is a demo decomposition of Rowan Smith's simulated top-down CO data image.

"""

from __future__ import division

import pickle
import os.path
import datetime
import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs
from astropy.io.fits import getdata, getheader
import astropy.io.fits as fits

from astrodendro.scatter import Scatter
from dendrogal.integrated_viewer import IntegratedViewer
from dendrogal.reid_distance_assigner import make_reid_distance_column
from dendrogal.assign_physical_values import assign_size_mass_alpha_pressure

data_path = os.path.expanduser("~/Dropbox/DendroGal/")

data, header = getdata(data_path+"WCO.fits", header=True)

d = astrodendro.Dendrogram.compute(data, verbose=True, min_value=1e-2, min_delta=1e-2, min_npix=100)

