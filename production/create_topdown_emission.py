"""
This is my first crack at making a "top down" EMISSION view of the Milky Way.

"""

from __future__ import division

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
