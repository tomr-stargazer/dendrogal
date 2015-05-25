"""
Turns the astropy.table catalog into a deluxetable for LaTeX.

"""

from __future__ import division
import os.path

import astropy.table

paper_directory = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def make_latex_table(input_table):

    filename = paper_directory+"test_table_thing.tex"

    latex_table = astropy.table.Table()

    # = input_table.copy(copy_data=True)

    latex_table['$l$ (\degr)'] = input_table['x_cen']
    latex_table['$b$ (\degr)'] = input_table['y_cen']
    latex_table['$v_{LSR}$ (km s$^{-1}$)'] = input_table['v_cen']

    latex_table['$\sigma_r$ (\degr)'] = input_table['radius']

    latex_table['Distance (kpc)'] = input_table['distance']
    latex_table['KDA resolution'] = input_table['KDA_resolution']

    latex_table['$R$ (pc)'] = input_table['size']
    latex_table['Mass ($M_\odot$)'] = input_table['mass']

    latex_table.write(filename, format='ascii.aastex')

