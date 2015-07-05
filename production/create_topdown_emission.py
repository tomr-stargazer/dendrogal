"""
This is my first crack at making a "top down" EMISSION view of the Milky Way.

Much of this is hackish and needs to be formalized into a class sometime, 
so that the various disjoint parts (e.g. "how big is the grid? 
what is the scaling between pc <-> px?") can be kept consistent.

Real (astropy) units need to be introduced too, as well as tests.

More reasons we should make it into a class: 
    * eventually we'll want to "downsample" the resulting image by 
      convolving it with different PSFs
    * it'd be nice to have a standard "now make a figure of it" method
      so we don't have to call plt.imshow all the time
    * Maybe we'll want to add noise or something to it

"""

from __future__ import division
import os.path

import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel

output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

# from production.dame_color_dict import dame_cmap

grid = np.zeros((2000, 2000))

def cloud_emission(mass, radius, center=(0,0)):

    this_clouds_emission = np.zeros_like(grid)

    # M = 4.4 X2 L_co
    L = mass / (4.4)

    F = L / radius**2

    # all the pixels within `radius` of `center`
    x_indices, y_indices = np.indices(grid.shape)

    x_offset = x_indices - center[0]
    y_offset = y_indices - center[1]

    cloud_circle = [(x_offset**2 + y_offset**2)**(1/2) < radius]

    this_clouds_emission[cloud_circle] += F

    return this_clouds_emission


# see definition here: http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
def gaussian_2d(x, y, A, sigma_r, x0, y0):
    r"""
    Two dimensional Gaussian with equal sigma_x, sigma_y.

    .. math:: f(x,y) = A \exp [ \frac{(x-x_0)^2+(y-y_0)^2}{2 \sigma_r^2} ]

    """

    return A * np.exp(-( (x-x0)**2 + (y-y0)**2)/(2*sigma_r**2) )


def cloud_emission_2d_gaussian(mass, radius_px, center_px=(0,0)):

    cloud_emission = np.zeros_like(grid)

    sigma_r = radius_px / 1.9

    A = mass / (4.4 * 2*np.pi * sigma_r**2)
    x0, y0 = center_px

    cloud_emission_function = lambda x, y: gaussian_2d(x, y, A, sigma_r, x0, y0)

    x_indices, y_indices = np.indices(cloud_emission.shape)

    cloud_emission += cloud_emission_function(x_indices, y_indices)

    return cloud_emission


def cheaply_convert_pc_to_px(radius_pc):

    # the grid is 25 kpc across.
    pc_per_px = 25000 / grid.shape[0] 

    radius_px = radius_pc / pc_per_px

    return radius_px


def transform_xgal_ygal_to_xgrid_ygrid(xgal, ygal):

    # the grid center is at xgal=0, ygal=0
    xgal0_in_px = grid.shape[0]/2
    ygal0_in_px = grid.shape[1]/2

    x_px = cheaply_convert_pc_to_px(xgal)
    x_grid = x_px + xgal0_in_px

    y_px = cheaply_convert_pc_to_px(ygal)
    y_grid = y_px + ygal0_in_px

    return x_grid, y_grid


def area_per_px():
    """ Returns the area of each pixel in pc^2. """

    len_per_px = cheaply_convert_pc_to_px(1)
    area = len_per_px**2

    return area


def len_per_px():
    """ Returns the area of each pixel in pc^2. """

    return cheaply_convert_pc_to_px(1)


def cloud_emission_from_catalog(catalog):

    beginning = datetime.datetime.now()

    cloud_emission_sum = np.zeros_like(grid)

    for i, row in enumerate(catalog):

        x_grid, y_grid = transform_xgal_ygal_to_xgrid_ygrid(row['x_gal']*1000, row['y_gal']*1000)
        radius_px = cheaply_convert_pc_to_px(row['size'])
        mass = row['mass']

        cloud_emission = cloud_emission_2d_gaussian(mass, radius_px, (x_grid, y_grid))
        if not np.isnan(cloud_emission).any():

            cloud_emission_sum += cloud_emission

        if i%10==0:
            print i
        # if i>100:
        #     break

    end = datetime.datetime.now()
    time_elapsed = (end - beginning)
    print " *** Image generation took {0}".format(time_elapsed)

    return cloud_emission_sum / area_per_px()


def image_convolution_and_display(catalog, emission=None, resolution_stddev=1.6):

    if emission is None:
        emission = cloud_emission_from_catalog(catalog)

    g_kernel = Gaussian2DKernel(stddev=resolution_stddev)

    result = convolve_fft(emission, g_kernel)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    image = ax.imshow(np.sqrt(result+1e-6), cmap='gist_heat_r', origin='lower', vmax=300)
    ax.invert_xaxis()
    ax.set_xticks([0*2, 250*2, 500*2, 750*2, 1000*2])
    ax.set_xticklabels([12.5, 6.25, 0, -6.25, -12.5])
    ax.set_yticks([0*2, 250*2, 500*2, 750*2, 1000*2])
    ax.set_yticklabels([-12.5, -6.25, 0, 6.25, 12.5])

    fig.colorbar(image)

    return fig


def make_emission_maps_for_paper(emission=None, catalog=None, save=False):

    if emission.shape != grid.shape:
        print 'size mismatch'
        raise ValueError("grid and emission are not same size - problems will occur")

    if emission is None:
        if catalog is None:

            from dendrogal.production.cloud_catalog_combiner import extract_and_combine_catalogs
            catalog = extract_and_combine_catalogs()

        emission = cloud_emission_from_catalog(catalog)

    paws_rez = cheaply_convert_pc_to_px(40)
    paws_10x_rez = cheaply_convert_pc_to_px(400)
    paws_5x_rez = cheaply_convert_pc_to_px(200)

    fig1 = image_convolution_and_display(None, emission, resolution_stddev=paws_rez)
    fig2 = image_convolution_and_display(None, emission, resolution_stddev=paws_10x_rez)
    fig3 = image_convolution_and_display(None, emission, resolution_stddev=paws_5x_rez)

    if save:
        fig1.savefig(output_path+"MW_sim_PAWS_40pc_rez.pdf", bbox_inches='tight')
        fig2.savefig(output_path+"MW_sim_400pc_rez.pdf", bbox_inches='tight')

    return fig1, fig2, fig3
