"""
Some python code to interface with Mark Reid's fortran tool.

Hopefully I can figure out how to do this in a non-boneheaded way.

"""

from __future__ import division

import subprocess
from subprocess import Popen, PIPE, STDOUT

import astropy.table

test_string = "lol    010203.04 121314.5 123 0"

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

astropy.table.Table.read(lines_of_data, format='ascii')

