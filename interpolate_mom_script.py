"""
A script that interpolates the moment-masked data and saves them to file.

"""

from __future__ import division

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
from astrodendro_analysis.integrated_viewer import IntegratedViewer
from astrodendro_analysis.reid_distance_assigner import make_reid_distance_column
from astrodendro_analysis.assign_physical_values import assign_size_mass_alpha_pressure

from dame_interpolation import interpolate_datacube

data_path = os.path.expanduser("~/Dropbox/College/Astro99/DATA/")

def create_intermom_file(filename, memmap=False, clobber=True):

    beginning = datetime.datetime.now()

    if "mom.fits" not in filename:
        raise ValueError("This function is only intended for files ending in *mom.fits")

    data, header = getdata(filename, memmap=memmap, header=True)

    new_filename = filename.rstrip("mom.fits") + "mominterp.fits"

    clobber_string = ""
    if os.path.isfile(new_filename) and clobber:
        clobber_string = "(clobbered)"
    elif os.path.isfile(new_filename):
        print "{0} not saved: clobber=False".format(new_filename)
        return

    new_data = interpolate_datacube(data)

    try:
        fits.writeto(new_filename, new_data, header, clobber=clobber)        
    except IOError, e:
        print "File not saved: {0}".format(e)
        return

    end = datetime.datetime.now()
    time_elapsed = (end - beginning)
    print "Time elapsed for {1}: {0} {2}".format(time_elapsed, new_filename.lstrip('DHT').rstrip("_mominterp.fits"), clobber_string)

    return


