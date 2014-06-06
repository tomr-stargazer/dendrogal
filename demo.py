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
from astrodendro_analysis.reid_distance_assigner import make_reid_distance_column
from astrodendro_analysis.assign_physical_values import assign_size_mass_alpha_pressure

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
    
def cogal_downsampled_demo(**kwargs):
    return downsampled_demo('COGAL_all_mom.fits', **kwargs)

def small_demo(**kwargs):
    return downsampled_demo('DHT17_Quad2_bw_mom.fits', **kwargs)

def orion_demo(**kwargs):
    return downsampled_demo('DHT27_Orion_mom.fits', **kwargs)

def downsampled_demo(data_file, downsample_factor=4, transpose_tuple=(2,0,1),
                     min_value=0.01, min_delta=0.005, min_npix=2000):

    df = downsample_factor
    tt = transpose_tuple

    print "loading data: ..."
    datacube, datacube_header = getdata(data_path+data_file, memmap=True,
                                        header=True)

    print "transposing, downsampling, and unit-converting data: ..."
    datacube_dt, datacube_dt_header = \
      downsample_and_transpose_data_and_header(datacube, datacube_header, df, tt)
    datacube_dt_wcs = wcs.wcs.WCS(datacube_dt_header)

    beam_size = 1/8 * u.deg

    # Convert the data from kelvin to jansky-per-beam
    omega_beam = np.pi * (0.5 * beam_size)**2 # Beam width is 1/8 degree
    frequency = 115 * u.GHz
    K_to_Jy = u.K.to(u.Jy, equivalencies=
                     u.brightness_temperature(omega_beam, frequency))
    datacube_dt_jansky_perbeam = datacube_dt * K_to_Jy 

    print "computing dendrogram: ..."
    d = astrodendro.Dendrogram.compute(
        datacube_dt_jansky_perbeam,
        min_value=min_value*K_to_Jy, min_delta=min_delta*K_to_Jy,  #these are arbitrary
        min_npix=min_npix//df**3, verbose=True)

    v_scale = datacube_dt_header['cdelt3']
    v_unit = u.km / u.s
    l_scale = datacube_dt_header['cdelt1']
    b_scale = datacube_dt_header['cdelt2']
    
    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam # According to A. Ginsburg
    metadata['spatial_scale'] = b_scale * u.deg
    metadata['velocity_scale'] = v_scale * v_unit
    metadata['wavelength'] = (c.c / frequency).to('mm')
    metadata['beam_major'] = beam_size
    metadata['beam_minor'] = beam_size    
    metadata['vaxis'] = 0 # keep it this way if you think the (post-downsample/transposed) input data is (l, b, v)
    metadata['wcs'] = datacube_dt_wcs

    catalog = astrodendro.ppv_catalog(d, metadata)

    if catalog['flux'].unit.is_equivalent('Jy'):
        # Workaround because flux is computed wrong

        flux = quantify_column(catalog['flux'])
        area_exact = catalog['area_exact'].unit*catalog['area_exact'].data

        flux_kelvin = flux.to('K', equivalencies=u.brightness_temperature(area_exact, frequency))

        flux_kelvin_kms_deg2 = flux_kelvin * metadata['velocity_scale'] * area_exact

        catalog.add_column(astropy.table.Column(data=flux_kelvin_kms_deg2, name='flux_kelvin_kms_deg2'))

    return d, catalog, datacube_dt_header, metadata

def multiple_linked_viewer_demo(demo=cogal_downsampled_demo, galactic=False, 
                                **kwargs):

    d, catalog, cogal_dt_header, metadata = demo(**kwargs)

    # DISTANCES
    try:
        reid = make_reid_distance_column(catalog)
        catalog['Distance'] = reid['D_k']        
    except Exception, e:
        print "DISTANCE ASSIGNMENT FAILED:", e
        catalog['Distance'] = 1000*np.ones_like(catalog['x_cen'])
        catalog['Distance'].unit = u.kpc

    # SIZE MASS VIRIAL
    s, m, v, p = assign_size_mass_alpha_pressure(catalog)

    catalog['size'] = astropy.table.Column(data=s, name='size')
    catalog['mass'] = astropy.table.Column(data=m, name='mass')
    catalog['virial'] = astropy.table.Column(data=v, name='virial')        
    catalog['pressure'] = astropy.table.Column(data=p, name='pressure')        

    if galactic:
        dv = d.viewer(galactic=True)
    else:
        dv = d.viewer()

    iv = IntegratedViewer(d, dv.hub)

    ds0 = Scatter(d, dv.hub, catalog, 'radius', 'v_rms')
    ds1 = Scatter(d, dv.hub, catalog, 'radius', 'flux')
    ds2 = Scatter(d, dv.hub, catalog, 'size', 'v_rms')
    ds3 = Scatter(d, dv.hub, catalog, 'mass', 'virial')
    ds4 = Scatter(d, dv.hub, catalog, 'size', 'pressure')

    scatter_viewers = [ds0, ds1, ds2, ds3, ds4]

    return_dict = {}
    return_dict['dendrogram'] = d
    return_dict['catalog'] = catalog
    return_dict['metadata'] = metadata
    return_dict['dendrogram_viewer'] = dv
    return_dict['integrated_viewer'] = iv
    return_dict['scatter_viewers'] = scatter_viewers

    return return_dict


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

def quantify_column(column):

    return (column.data * column.unit)