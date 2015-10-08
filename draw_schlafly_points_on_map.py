"""
Our goal here is to draw the Schlafly points on a map.

"""

from __future__ import division
import os

import numpy as np
from astropy.table import Table
import astropy.wcs as wcs
# load up the Schlafly data

from dendrogal.overlay_dame_cartoon import corners, px_per_degree

code_directory = os.path.dirname(__file__)
schlafly_path = os.path.join(code_directory, "schlafly")

big_cloud_table = Table.read( os.path.join(schlafly_path, "schafly_bigcloud_table.fit") )
mbm_cloud_table = Table.read( os.path.join(schlafly_path, "schafly_mbmcloud_table.fit") )

def draw_cloud_table_on_wcsaxes(ax, cloud_table):

    wcs_object= ax.wcs.sub([wcs.WCSSUB_CELESTIAL])
    wcs_object.wcs.bounds_check(False,False)

    l_coords = cloud_table['GLON']
    b_coords = cloud_table['GLAT']

    output = wcs_object.wcs_world2pix(np.vstack((l_coords, b_coords)).T,0)

    l_px = output[:,0]
    b_px = output[:,1]

    lines = ax.plot(l_px, b_px, 'm.', ms=8)

    return lines


def draw_big_cloud_table(ax):
    return draw_cloud_table_on_wcsaxes(ax, big_cloud_table)


def draw_mbm_cloud_table(ax):
    return draw_cloud_table_on_wcsaxes(ax, mbm_cloud_table)