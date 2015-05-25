"""
Turns the astropy.table catalog into a deluxetable for LaTeX.

"""

from __future__ import division

import os.path

paper_directory = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def make_latex_table(input_table):

    filename = paper_directory+"test_table_thing.tex"

    input_table.write(filename, format='ascii.aastex')

