"""
A script that interpolates the moment-masked data and saves them to file.

"""

from __future__ import division

import os.path
import datetime
import numpy as np

from astropy import wcs
from astropy.io.fits import getdata, getheader
import astropy.io.fits as fits

from dame_interpolation import interpolate_datacube
from dame_moment_masking import moment_mask

data_path = os.path.expanduser("~/Dropbox/College/Astro99/DATA/")

def create_intermom_file(filename, memmap=False, clobber=False, rms_noise=None):

    beginning = datetime.datetime.now()

    print "Beginning {0} at {1}".format(filename, datetime.datetime.strftime(beginning,"%Y-%m-%d %H:%M:%S"))

    if "_mom.fits" in filename:

        data, header = getdata(filename, memmap=memmap, header=True)

        new_filename = filename.rstrip("mom.fits") + "mominterp.fits"

        clobber_string = ""
        if os.path.isfile(new_filename) and clobber:
            clobber_string = "(clobbered)"
        elif os.path.isfile(new_filename):
            print "{0} not saved: clobber=False".format(new_filename)
            return

        new_data = interpolate_datacube(data)

    elif "_interp.fits" in filename: 
        data, header = getdata(filename, memmap=memmap, header=True)

        if rms_noise is None:
            raise ValueError("Please provide an rms_noise.")
        new_filename = filename.rstrip("interp.fits") + "interp_mom.fits"

        clobber_string = ""
        if os.path.isfile(new_filename) and clobber:
            clobber_string = "(clobbered)"
        elif os.path.isfile(new_filename):
            print "{0} not saved: clobber=False".format(new_filename)
            return

        new_data = moment_mask(data, rms_noise)

    else:
        raise ValueError("This function is only intended for files ending in *mom.fits or *interp.fits")

    try:
        fits.writeto(new_filename, new_data, header, clobber=clobber)        
    except IOError, e:
        print "File not saved: {0}".format(e)
        return

    end = datetime.datetime.now()
    time_elapsed = (end - beginning)
    print "Time elapsed for {1}: {0} {2}".format(time_elapsed, new_filename.lstrip('DHT').replace("_mominterp.fits", ""), clobber_string)

    return

def create_intermom_folder(path=data_path, **kwargs):

    file_list = os.listdir(path)

    mom_files = [x for x in file_list if 'mom.fits' in x]

    for filename in mom_files:

        create_intermom_file(path+filename, **kwargs)
