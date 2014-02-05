"""
Some python code to interface with Mark Reid's fortran tool.

Hopefully I can figure out how to do this in a non-boneheaded way.

"""

from __future__ import division

import subprocess
from subprocess import Popen, PIPE, STDOUT

import astropy.table

test_string = "test    010203.04 121314.5 -10 1"

f = open('source_file.dat', 'w')
f.write(test_string)
f.close()

p = Popen(
    ['/Users/tsrice/Documents/Code/mark_reid_kdist/revised_kinematic_distance'],
    shell = False, stdin=PIPE, stdout=PIPE)

output, error = p.communicate()

print output
print error

lines_of_data = [x for x in output.split('\n') if '!' not in x]

kd_output = astropy.table.Table.read(lines_of_data, format='ascii')

"""
! Source     Gal Long  Gal Lat    V_lsr     V_rev    Rev. D_k     +/-
!              (deg)    (deg)    (km/s)    (km/s)     (kpc)      (kpc)
"""

kd_output.rename_column('col1', 'Source')
kd_output.rename_column('col2', 'gal_long')
kd_output.rename_column('col3', 'gal_lat')
kd_output.rename_column('col4', 'V_lsr')
kd_output.rename_column('col5', 'V_rev')
kd_output.rename_column('col6', 'D_k')
kd_output.rename_column('col7', 'error_D_k_plus')
kd_output.rename_column('col8', 'error_D_k_minus')

print kd_output
