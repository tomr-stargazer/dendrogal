"""
This is a demo to open, downsample, and dendrogram some sample data.

The data are allowed to be the "actual" thing we're trying to dendrogram I guess.

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
try:
    from dendrogal.reid_distance_assigner import make_reid_distance_column
    from dendrogal.assign_physical_values import assign_size_mass_alpha_pressure
except:
    pass

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

def resample_3d_variable(array, x_factor, y_factor, z_factor):

    nx, ny, nz = np.shape(array)

    nx_new = nx // x_factor
    ny_new = ny // y_factor
    nz_new = nz // z_factor    

    array2 = np.zeros((nx_new, ny, nz))
    for i in range(nx_new):
        array2[i, :, :] = np.nanmean(array[i * x_factor:(i + 1) * x_factor, :, :], axis=0)
    del array

    array3 = np.zeros((nx_new, ny_new, nz))
    for j in range(ny_new):
        array3[:, j, :] = np.nanmean(array2[:, j * y_factor:(j + 1) * y_factor, :], axis=1)
    del array2

    array4 = np.zeros((nx_new, ny_new, nz_new))
    for k in range(nz_new):
        array4[:, :, k] = np.nanmean(array3[:, :, k * z_factor:(k + 1) * z_factor], axis=2)
    del array3

    return array4

def recenter_wcs_header(input_header, central_value=0):
    """ 
    Sets the header CRVAL on zero if it's not already. 

    A quick-and-dirty implementation.

    """

    new_header = input_header.copy()

    # in principle we could iterate through dimensions but.... not yet
    if input_header['CRVAL1'] != central_value:
        # how do we figure out where it hits zero? Use CDELT1 and do some math?
        degrees_correction = central_value - input_header['CRVAL1']
        pixels_correction = degrees_correction / input_header['CDELT1']

        new_header['CRVAL1'] = central_value
        new_header['CRPIX1'] = input_header['CRPIX1'] + pixels_correction

    return new_header


def downsample_and_transpose_data_and_header(input_data, input_header, 
                                             downsample_factor=4, 
                                             transpose_tuple=(2,0,1), resample=False, recenter=True):
    
    df = downsample_factor
    tt = transpose_tuple

    # this bit is to support resample_3d_variable
    if type(df) is not int:
        if len(df) != len(input_data.shape):
            raise ValueError("downsample_factor invalid")

        new_data = resample_3d_variable(input_data, *df).transpose(*tt)
        df = tuple(df[i] for i in tt)

    elif resample and df > 1:
        new_data = resample_3d(input_data, df).transpose(*tt)
        df = (df, df, df)
    else:
        # Someday I may improve this with something less crude.
        new_data = input_data[::df, ::df, ::df].transpose(*tt)
        df = (df, df, df)

    new_header = input_header.copy()

    # let's transpose. and downsample.
    new_header['naxis'+str(tt[0]+1)] = input_header['naxis1'] // df[0]
    new_header['naxis'+str(tt[1]+1)] = input_header['naxis2'] // df[1]
    new_header['naxis'+str(tt[2]+1)] = input_header['naxis3'] // df[2] 

    new_header['ctype'+str(tt[0]+1)] = input_header['ctype1']
    new_header['ctype'+str(tt[1]+1)] = input_header['ctype2']
    new_header['ctype'+str(tt[2]+1)] = input_header['ctype3']

    new_header['crval'+str(tt[0]+1)] = input_header['crval1']
    new_header['crval'+str(tt[1]+1)] = input_header['crval2']
    new_header['crval'+str(tt[2]+1)] = input_header['crval3']

    new_header['cdelt'+str(tt[0]+1)] = input_header['cdelt1'] * df[0]
    new_header['cdelt'+str(tt[1]+1)] = input_header['cdelt2'] * df[1]
    new_header['cdelt'+str(tt[2]+1)] = input_header['cdelt3'] * df[2]

    if resample:
        new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df[0]
        new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df[1]
        new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df[2]
    else:
        new_header['crpix'+str(tt[0]+1)] = (input_header['crpix1'] - 1)//df[0] + 1
        new_header['crpix'+str(tt[1]+1)] = (input_header['crpix2'] - 1)//df[1] + 1
        new_header['crpix'+str(tt[2]+1)] = (input_header['crpix3'] - 1)//df[2] + 1

    if recenter:
        return_header = recenter_wcs_header(new_header)
    else:
        return_header = new_header

    return new_data, return_header

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

# After trial and error I found that axes=2 is the longitude axis when it comes to periodic_neighbors, despite what the docs say
longitude_neighbors = astrodendro.dendrogram.periodic_neighbours(2)

def cogal_downsampled_demo(**kwargs):
    return downsampled_demo('COGAL_all_mom.fits', 
                            neighbours=longitude_neighbors,
                            **kwargs)

def cogal_resampled2_demo(downsample_factor=1, transpose_tuple=(0,1,2), **kwargs):
    return downsampled_demo('COGAL_all_mom_downsampled_by_2.fits', transpose_tuple=transpose_tuple, 
                            downsample_factor=downsample_factor, neighbours=longitude_neighbors, **kwargs)

def cogal_resampled_demo(resample, downsample_factor=1, transpose_tuple=(0,1,2), min_npix=None, **kwargs):
    return downsampled_demo('COGAL_all_mom_downsampled_by_{0}.fits'.format(resample), 
                            transpose_tuple=transpose_tuple, downsample_factor=downsample_factor,
                            neighbours=longitude_neighbors, min_npix=min_npix or 2000//resample, **kwargs)

def cogal_local_demo(**kwargs):
    return downsampled_demo('COGAL_local_mom.fits', downsample_factor=2, neighbours=longitude_neighbors)

def cogal_local_resampled_demo(resample, downsample_factor=1, transpose_tuple=(0,1,2), min_npix=None, **kwargs):
    return downsampled_demo('COGAL_local_mom_downsampled_by_{0}.fits'.format(resample), 
                            transpose_tuple=transpose_tuple, downsample_factor=downsample_factor,
                            neighbours=longitude_neighbors, min_npix=min_npix or 2000//resample, **kwargs)

def cogal_deep_demo(downsample_factor=2, mom_interp='mom', **kwargs):
    return downsampled_demo('COGAL_deep_{0}.fits'.format(mom_interp), downsample_factor=downsample_factor, neighbours=longitude_neighbors, **kwargs)

def cogal_deep_resampled_demo(resample, downsample_factor=1, transpose_tuple=(0,1,2), min_npix=None, **kwargs):
    return downsampled_demo('COGAL_deep_mom_downsampled_by_{0}.fits'.format(resample), 
                            transpose_tuple=transpose_tuple, downsample_factor=downsample_factor,
                            neighbours=longitude_neighbors, min_npix=min_npix or 2000//resample, **kwargs)

def small_demo(**kwargs):
    return downsampled_demo('DHT17_Quad2_bw_mom.fits', **kwargs)

def orion_demo(**kwargs):
    return downsampled_demo('DHT27_Orion_mom.fits', **kwargs)

def orion_demo_mominterp(**kwargs):
    return downsampled_demo('DHT27_Orion_mominterp.fits', **kwargs)

def chamaeleon_demo(**kwargs):
    return downsampled_demo('DHT34_Chamaeleon_mom.fits', **kwargs)

def perseus_demo(**kwargs):
    return downsampled_demo('DHT21_Taurus_mom.fits', **kwargs)

def polaris_demo(**kwargs):
    return downsampled_demo("DHT16_Polaris_bw_mom.fits", **kwargs)

def ophiuchus_demo(**kwargs):
    return downsampled_demo("DHT37_Ophiuchus_mom.fits", **kwargs)

def firstquad_demo(**kwargs):
    return downsampled_demo("DHT08_Quad1_mom.fits", **kwargs)

def firstquad_demo_mominterp(**kwargs):
    return downsampled_demo("DHT08_Quad1_mominterp.fits", **kwargs)

def secondquad_demo(**kwargs):
    return downsampled_demo("DHT17_Quad2_bw_mom.fits", **kwargs)

def secondquad_demo_mominterp(**kwargs):
    return downsampled_demo("DHT17_Quad2_bw_mominterp.fits", **kwargs)

def thirdquad_demo_mominterp(**kwargs):
    return downsampled_demo("DHT31_Quad3_mominterp.fits", **kwargs)


def downsampled_demo(data_file, downsample_factor=4, transpose_tuple=(2,0,1),
                     min_value=0.01, min_delta=0.005, min_npix=2000, 
                     neighbours=None, resample=False, compute_catalog=True, recenter=True,
                     data_path=data_path, memmap=True):

    df = downsample_factor
    tt = transpose_tuple

    print("\n ** Downsampled Demo Parameters:"
          "\n    downsample_factor={0}"
          "\n    transpose_tuple={1}"
          "\n    min_value={2}"
          "\n    min_delta={3}"
          "\n    min_npix={4}\n".format(downsample_factor, transpose_tuple, min_value, min_delta, min_npix))

    beginning = datetime.datetime.now()
    beginning_string = datetime.datetime.strftime(beginning,"%Y-%m-%d %H:%M:%S")

    print("\n ** Beginning at: {0}".format(beginning_string))

    print "\n** loading data from {0}: ...\n".format(data_path+data_file)
    datacube, datacube_header = getdata(data_path+data_file, memmap=memmap,
                                        header=True)

    print "\n** transposing {0} and downsampling ({1}) data: ...\n".format(tt, df)
    datacube_dt, datacube_dt_header = \
      downsample_and_transpose_data_and_header(datacube, datacube_header, df, tt, resample=resample, recenter=recenter)
    datacube_dt_wcs = wcs.wcs.WCS(datacube_dt_header)

    datacube_dt_wcs.wcs.bounds_check(pix2world=False)

    beam_size = 1/8 * u.deg
    frequency = 115 * u.GHz

    print "\n** computing dendrogram on data with dimensions {} and n_pixels {:,}: ...\n".format(np.shape(datacube_dt), np.size(datacube_dt))
    d = astrodendro.Dendrogram.compute(
        datacube_dt,
        min_value=min_value, min_delta=min_delta,  #these are arbitrary
        min_npix=min_npix, verbose=True, neighbours=neighbours, wcs=datacube_dt_wcs)

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

    print "\n** computing catalog for {:,} structures: ...\n".format(len(d))
    if compute_catalog:
        catalog = astrodendro.ppv_catalog(d, metadata, verbose=True)
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

    end = datetime.datetime.now()
    end_string = datetime.datetime.strftime(end,"%Y-%m-%d %H:%M:%S")

    print("\n ** Ending at: {0}".format(end_string))

    time_elapsed = (end - beginning)
    print "\n ** Time elapsed: {0}".format(time_elapsed)

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


savepath = os.path.expanduser("~/Documents/Code/dendrogal/saved_dendrogram/")
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

def resample_transpose_and_save(data_file, downsample_factor=2, transpose_tuple=(2,0,1), output_path=None, clobber=True, **kwargs):

    df = downsample_factor
    tt = transpose_tuple

    print "loading data: ..."
    datacube, datacube_header = getdata(data_path+data_file, memmap=True,
                                        header=True)

    print "SHAPE of input data: {0}".format(np.shape(datacube))

    print "transposing and resampling data: ..."
    datacube_dt, datacube_dt_header = \
      downsample_and_transpose_data_and_header(datacube, datacube_header, df, tt, resample=True)

    print "SHAPE of output data: {0}".format(np.shape(datacube_dt))
    expected_shape = [x/df for x in np.shape(datacube)]
    expected_shape = expected_shape[tt[0]], expected_shape[tt[1]], expected_shape[tt[2]]
    print "EXPECTED SHAPE of output data: {0}".format(expected_shape)

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
            return datacube, datacube_header, datacube_dt, datacube_dt_header

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

def reduce_selection_to_principal_branches(structures):
    """
    Remove all degenerate structures in a list 

    (i.e., structures that are descendants of other structures) 

    """

    structure_set = set(structures)

    # sort by size
    structures_by_descendants = sorted(structures, key=lambda k: len(k.descendants))

    # from smallest to largest:
    for little_struct, i in zip(structures_by_descendants, range(len(structures_by_descendants))):

        # ask from largest to smallest:
        for big_struct in structures_by_descendants[:i:-1]: # backwards, and skip everyone smaller than little_struct

            # "are you my ancestor?" (# wait, why isn't there an "ancestors" property?)
            if little_struct in big_struct.descendants: 

                # if so, remove the little guy
                structure_set.remove(little_struct)
                # and stop asking
                break 

    return list(structure_set)