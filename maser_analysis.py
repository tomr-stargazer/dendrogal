"""
Currently a mockup module just to see if we can compare masers to our dendrogram structures.

"""

import os

from astropy.coordinates import SkyCoord
from astropy.table import Table

table_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/tables/")

maser_table = Table.read(table_path+'Reid2014_parallaxes_apj490685t1_mrt.txt', format='ascii.cds')

signed_degrees = maser_table['DEd'].copy()
signed_degrees[maser_table['DE-'] == '-'] *= -1

ra_tuple = (maser_table['RAh'], maser_table['RAm'], maser_table['RAs'])
dec_tuple = (signed_degrees, maser_table['DEm'], maser_table['DEs'])

maser_coordinates = SkyCoord(ra=ra_tuple, dec=dec_tuple, unit=('hour', 'deg'))

# Code to figure out what structure a maser thing is on!

# Find pixel co-ordinates of click

def put_masers_in_dendrogram(dendrogram, catalog, header, metadata):

	y_wcs =  metadata['wcs']

	for coordinate, velocity in zip(maser_coordinates.galactic, maser_table['VLSR']):

		l_px, b_px, v_px = y_wcs.wcs_world2pix(coordinate.l, coordinate.b, velocity, 1)

		indices = (v_px, b_px, l_px) # indexing is funny

		try:
			structure = dendrogram.structure_at(indices)

			print structure

		except Exception, e:

			print e



# indices = (self.slice, iy, ix)

# # Select the structure
# structure = 
