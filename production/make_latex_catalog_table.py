"""
Turns the astropy.table catalog into a deluxetable for LaTeX.

"""

from __future__ import division
import os.path

import numpy as np

import astropy.table
import astropy.units as u

paper_directory = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")

def make_latex_table(input_table, start=0, end=27):

    input_table.sort('x_cen')

    filename = paper_directory+"test_table_thing.tex"

    latex_table = astropy.table.Table()

    # = input_table.copy(copy_data=True)

    caption_string_text = """ 
    Table of the {0} clouds identified in this work. 
    Column "KDA" refers to the resolution of the kinematic distance ambiguity; 
    'U' refers to unambiguous distances, while 'N' and 'F' denote near and far.
    The full table appears in the online version of the journal.
    """.format(len(input_table))

    caption_string = "\caption{"+caption_string_text+"}"

    latexdict = {}
    latexdict['caption'] = "Cloud Catalog"
    latexdict['tablefoot'] = caption_string

    latex_table['$l$'] = input_table['x_cen']
    latex_table['$l$'].unit = u.deg
    latex_table['$l$'].format = '%.2f'
    latex_table['$b$'] = input_table['y_cen']
    latex_table['$b$'].unit = u.deg
    latex_table['$b$'].format = '%.2f'
    latex_table[r'$v_{\textrm{LSR}}$'] = input_table['v_cen']
    latex_table[r'$v_{\textrm{LSR}}$'].unit = u.km * u.s**(-1)
    latex_table[r'$v_{\textrm{LSR}}$'].format = '%.2f'

    latex_table['$\sigma_r$'] = input_table['radius']
    latex_table['$\sigma_r$'].format = '%.2f'

    latex_table['$\sigma_v$'] = input_table['v_rms']
    latex_table['$\sigma_v$'].format = '%.2f'    

    latex_table['Distance'] = input_table['distance']
    latex_table['Distance'].format = '%.2f'    
    latex_table['KDA'] = input_table['KDA_resolution']

    latex_table['$R$'] = input_table['size']
    latex_table['$R$'].format = '%.2f'

    mass_string_list = [r"${0:.2f} \times 10^{1}$".format(mass/10**(np.floor(np.log10(mass))), int(np.floor(np.log10(mass)))) for mass in input_table['mass']]

    latex_table['Mass'] = mass_string_list
    latex_table['Mass'].unit = u.Msun

    latex_table[start:end].write(filename, format='ascii.aastex', latexdict=latexdict)

