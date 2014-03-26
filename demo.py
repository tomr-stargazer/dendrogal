"""
This is a demo to open, downsample, and dendrogram some sample data.

The data are allowed to be the "actual" thing we're trying to dendrogram I guess.

"""

from __future__ import division

import pickle
import numpy as np

import astropy
import astrodendro
import astropy.units as u
import astropy.constants as c

from astropy import wcs
from astropy.io.fits import getdata, getheader

from astrodendro.scatter import Scatter
from astrodendro_analysis.integrated_viewer import IntegratedViewer

data_path = "/Users/tsrice/Dropbox/college/Astro99/DATA/"

def downsample_and_transpose_data_and_header(input_data, input_header, 
                                             downsample_factor=4, 
                                             transpose_tuple=(2,0,1)):
    
    df = downsample_factor
    tt = transpose_tuple

    # Someday I may improve this with something less crude.
    new_data = input_data[::df, ::df, ::df].transpose(*tt)

    new_header = input_header.copy()

    # let's transpose. and downsample.
    new_header['naxis'+str(tt[0]+1)] = input_header['naxis1'] // df
    new_header['naxis'+str(tt[1]+1)] = input_header['naxis2'] // df
    new_header['naxis'+str(tt[2]+1)] = input_header['naxis3'] // df    

    new_header['ctype'+str(tt[0]+1)] = input_header['ctype1']
    new_header['ctype'+str(tt[1]+1)] = input_header['ctype2']
    new_header['ctype'+str(tt[2]+1)] = input_header['ctype3']

    new_header['crval'+str(tt[0]+1)] = input_header['crval1']
    new_header['crval'+str(tt[1]+1)] = input_header['crval2']
    new_header['crval'+str(tt[2]+1)] = input_header['crval3']

    new_header['cdelt'+str(tt[0]+1)] = input_header['cdelt1'] * df
    new_header['cdelt'+str(tt[1]+1)] = input_header['cdelt2'] * df
    new_header['cdelt'+str(tt[2]+1)] = input_header['cdelt3'] * df

    new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df + 1
    new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df + 1
    new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df + 1

    return new_data, new_header
    

def cogal_downsampled_demo(downsample_factor=4, transpose_tuple=(2,0,1)):

    df = downsample_factor
    tt = transpose_tuple

    print "loading data: ..."
    cogal, cogal_header = getdata(data_path+'COGAL_all_mom.fits', memmap=True,
                                  header=True)

    print "transposing, downsampling, and unit-converting data: ..."
    cogal_dt, cogal_dt_header = \
      downsample_and_transpose_data_and_header(cogal, cogal_header, df, tt)      
    cogal_dt_wcs = wcs.wcs.WCS(cogal_dt_header)

    beam_size = 1/8 * u.deg

    # Convert the data from kelvin to jansky-per-beam
    omega_beam = np.pi * (0.5 * beam_size)**2 # Beam width is 1/8 degree
    frequency = 115 * u.GHz
    K_to_Jy = u.K.to(u.Jy, equivalencies=
                     u.brightness_temperature(omega_beam, frequency))
    cogal_dt_jansky_perbeam = cogal_dt * K_to_Jy 

    print "computing dendrogram: ..."
    d = astrodendro.Dendrogram.compute(
        cogal_dt_jansky_perbeam,
        min_value=0.01*K_to_Jy, min_delta=0.005*K_to_Jy,  #these are arbitrary
        min_npix=2000//df**3, verbose=True)

    v_scale = cogal_dt_header['cdelt3']
    v_unit = u.km / u.s
    l_scale = cogal_dt_header['cdelt1']
    b_scale = cogal_dt_header['cdelt2']
    
    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam # According to A. Ginsburg
    metadata['spatial_scale'] = b_scale * u.deg
    metadata['velocity_scale'] = v_scale * v_unit
    metadata['wavelength'] = (c.c / frequency).to('mm')
    metadata['beam_major'] = beam_size
    metadata['beam_minor'] = beam_size    
    metadata['vaxis'] = 0 # keep it this way if you think the (post-downsample/transposed) input data is (l, b, v)
    metadata['wcs'] = cogal_dt_wcs

    catalog = astrodendro.ppv_catalog(d, metadata)

    return d, catalog, cogal_dt_header, metadata

def multiple_linked_viewer_demo(**kwargs):

    d, catalog, cogal_dt_header, metadata = cogal_downsampled_demo(**kwargs)

    dv = d.viewer()

    iv = IntegratedViewer(d, dv.hub)

    ds = Scatter(d, dv.hub, catalog, 'radius', 'flux')

    """ # this don't work -- the windows don't rescale when the figures do
    # just to get things to tile nicely on my screen...
    dv.fig.set_size_inches(16, 4)
    dv.fig.canvas.draw()
    ds.fig.set_size_inches(6, 4.4)
    ds.fig.canvas.draw()    """

    return dv, iv, ds


savepath = "/Users/tsrice/Documents/Code/astrodendro_analysis/saved_dendrogram/"
dendro_fame = "saved_dendrogram_object.hdf5"
catalog_fname = "saved_catalog_table.fits"
header_fname = "saved_header.fits"
metadata_fname = "saved_metadata.p"

def write_cogal_demo_to_file(**kwargs):

    d, catalog, cogal_dt_header, metadata = cogal_downsampled_demo(**kwargs)

    d.save_to(savepath+dendro_fame)
    catalog.write(savepath+catalog_fname)
    cogal_dt_header.tofile(savepath+header_fname)
    pickle.dump(metadata, open(savepath+metadata_fname, 'wb'))

    return None

# I noticed that "loaded" dendrograms tend not to respond well to click events... 
# the branches seem to lose their self._dendrogram attachments and everything breaks.
def load_demo_from_file(override=False):

    if not override:
        raise Exception("Don't use this method! The dendrograms fail to load properly.")

    d = astrodendro.Dendrogram.load_from(savepath+dendro_fame)
    catalog = astropy.table.Table.read(savepath+catalog_fname)
    cogal_dt_header = getheader(savepath+header_fname)
    metadata = pickle.load(open(savepath+metadata_fname, 'rb'))

    return d, catalog, cogal_dt_header, metadata