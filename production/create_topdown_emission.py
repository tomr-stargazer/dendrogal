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

import datetime

import numpy as np

# from production.dame_color_dict import dame_cmap

grid = np.zeros((1000, 1000))

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

    # the grid is 30 kpc across.
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

    return cloud_emission_sum