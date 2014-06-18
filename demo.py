"""
This is a demo to open, downsample, and dendrogram some sample data.

The data are allowed to be the "actual" thing we're trying to dendrogram I guess.

"""

from __future__ import division

import pickle
import os.path
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

data_path = os.path.expanduser("~/Dropbox/College/Astro99/DATA/")

# Downsamples a datacube by averaging. Code derived from aplpy.image_util (written by Tom Robitaille)
def resample_3d(array, factor):

    nx, ny, nz = np.shape(array)

    nx_new = nx // factor
    ny_new = ny // factor
    nz_new = nz // factor    

    array2 = np.zeros((nx_new, ny, nz))
    for i in range(nx_new):
        array2[i, :, :] = np.nanmean(array[i * factor:(i + 1) * factor, :, :], axis=0)
    del array

    array3 = np.zeros((nx_new, ny_new, nz))
    for j in range(ny_new):
        array3[:, j, :] = np.nanmean(array2[:, j * factor:(j + 1) * factor, :], axis=1)
    del array2

    array4 = np.zeros((nx_new, ny_new, nz_new))
    for k in range(nz_new):
        array4[:, :, k] = np.nanmean(array3[:, :, k * factor:(k + 1) * factor], axis=2)
    del array3

    return array4


def downsample_and_transpose_data_and_header(input_data, input_header, 
                                             downsample_factor=4, 
                                             transpose_tuple=(2,0,1), resample=False):
    
    df = downsample_factor
    tt = transpose_tuple

    if resample and df > 1:
        new_data = resample_3d(input_data, df).transpose(*tt)
    else:
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

    if resample:
        new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df
        new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df
        new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df
    else:
        new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df + 1
        new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df + 1
        new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df + 1

    return new_data, new_header

def make_2d_wcs_from_3d_wcs(input_wcs):
    """ Turns a PPV WCS object into a PP WCS object by modifying its header representation. """

    header = input_wcs.to_header()
    if header['WCSAXES'] != 3:
        raise ValueError("This function is only intended to work when WCSAXES == 3.")

    header_keys = header.keys()

    # strip the third axes everywhere
    for key in header_keys:
        if '3' in key:
            del header[key]

    # then declare there are only 2 axes
    header['WCSAXES'] = 2

    new_wcs = wcs.wcs.WCS(header)

    return new_wcs
    
def cogal_downsampled_demo(**kwargs):
    return downsampled_demo('COGAL_all_mom.fits', **kwargs)

def cogal_resampled2_demo(downsample_factor=1, **kwargs):
    return downsampled_demo('COGAL_all_mom_downsampled_by_2.fits', downsample_factor=downsample_factor, **kwargs)

def small_demo(**kwargs):
    return downsampled_demo('DHT17_Quad2_bw_mom.fits', **kwargs)

def orion_demo(**kwargs):
    return downsampled_demo('DHT27_Orion_mom.fits', **kwargs)

def perseus_demo(**kwargs):
    return downsampled_demo('DHT21_Taurus_mom.fits', **kwargs)

def ophiuchus_demo(**kwargs):
    return downsampled_demo('DHT37_Ophiuchus_mom.fits', **kwargs)

def downsampled_demo(data_file, downsample_factor=4, transpose_tuple=(2,0,1),
                     min_value=0.01, min_delta=0.005, min_npix=2000, resample=False, compute_catalog=True):

    df = downsample_factor
    tt = transpose_tuple

    print "loading data: ..."
    datacube, datacube_header = getdata(data_path+data_file, memmap=True,
                                        header=True)

    print "transposing, downsampling, and unit-converting data: ..."
    datacube_dt, datacube_dt_header = \
      downsample_and_transpose_data_and_header(datacube, datacube_header, df, tt, resample=resample)
    datacube_dt_wcs = wcs.wcs.WCS(datacube_dt_header)

    beam_size = 1/8 * u.deg
    frequency = 115 * u.GHz

    print "computing dendrogram: ..."
    d = astrodendro.Dendrogram.compute(
        datacube_dt,
        min_value=min_value, min_delta=min_delta,  #these are arbitrary
        min_npix=min_npix//df**3, verbose=True)

    v_scale = datacube_dt_header['cdelt3']
    v_unit = u.km / u.s
    l_scale = datacube_dt_header['cdelt1']
    b_scale = datacube_dt_header['cdelt2']

    metadata = {}
    metadata['data_unit'] = u.K
    metadata['spatial_scale'] = b_scale * u.deg
    metadata['velocity_scale'] = v_scale * v_unit
    metadata['wavelength'] = frequency # formerly: (c.c / frequency).to('mm') but now compute_flux can handle frequency in spectral equivalency
    metadata['beam_major'] = beam_size
    metadata['beam_minor'] = beam_size    
    metadata['vaxis'] = 0 # keep it this way if you think the (post-downsample/transposed) input data is (l, b, v)
    metadata['wcs'] = datacube_dt_wcs

    print "computing catalog: ..."
    if compute_catalog:
        catalog = astrodendro.ppv_catalog(d, metadata)
    else:
        catalog = None
        return d, catalog, datacube_dt_header, metadata

    if catalog['flux'].unit.is_equivalent('Jy'):
        # Workaround because flux is computed wrong

        flux = quantify_column(catalog['flux'])
        area_exact = catalog['area_exact'].unit*catalog['area_exact'].data

        # average brightness temperature integrated over area_exact
        flux_kelvin = flux.to('K', equivalencies=u.brightness_temperature(area_exact, frequency))
        # flux integrated over area and velocity
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


savepath = os.path.expanduser("~/Documents/Code/astrodendro_analysis/saved_dendrogram/")
dendro_fame = "saved_dendrogram_object.hdf5"
catalog_fname = "saved_catalog_table.fits"
header_fname = "saved_header.fits"
metadata_fname = "saved_metadata.p"

def write_cogal_demo_to_file(**kwargs):

    if (os.path.isfile(savepath+dendro_fame) or os.path.isfile(savepath+catalog_fname) or
        os.path.isfile(savepath+header_fname) or os.path.isfile(savepath+metadata_fname)):

        raise Exception("Files already exist! we should delete them.")

    d, catalog, cogal_dt_header, metadata = cogal_downsampled_demo(**kwargs)

    d.save_to(savepath+dendro_fame)
    catalog.write(savepath+catalog_fname)
    cogal_dt_header.tofile(savepath+header_fname)
    pickle.dump(metadata, open(savepath+metadata_fname, 'wb'))

    return None

def resample_transpose_and_save(data_file, downsample_factor=2, transpose_tuple=(2,0,1), output_path=None, clobber=False, **kwargs):

    df = downsample_factor
    tt = transpose_tuple

    print "loading data: ..."
    datacube, datacube_header = getdata(data_path+data_file, memmap=True,
                                        header=True)

    print "transposing, downsampling, and unit-converting data: ..."
    datacube_dt, datacube_dt_header = \
      downsample_and_transpose_data_and_header(datacube, datacube_header, df, tt, resample=True)

    print "saving data: ..."
    output_path = output_path or "{0}{1}_downsampled_by_{2}.fits".format(data_path, data_file.rstrip('.fits'), df)

    try:
        fits.writeto(output_path, datacube_dt, datacube_dt_header)
    except IOError, e:
        if clobber:
            os.remove(output_path)
            fits.writeto(output_path, datacube_dt, datacube_dt_header)
        else:
            print "File not saved: {0}".format(e)
            return

    print "Resampled file saved to {0}".format(output_path)

    return output_path

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