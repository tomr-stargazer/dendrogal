"""
Some python code to interface with Mark Reid's fortran tool.

Hopefully I can figure out how to do this in a non-boneheaded way.

"""

from __future__ import division

import subprocess

from subprocess import Popen, PIPE, STDOUT

#p = Popen(['grep', 'f'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)

#grep_stdout = p.communicate(input='one\ntwo\nthree\nfour\nfive\nsix\n')[0]

#StdinCommand = '''
#    MolPro code
#'''

# given some kind of input variables, generate the string or file 
# that revised_kinematic_distance.f needs

#   name rastring destring vstring fstring
test_string = "lol    010203.04 121314.5 123 0"

p = Popen(['/Users/tsrice/Documents/Code/mark_reid_kdist/revised_kinematic_distance'], 
          shell = False, stdin = PIPE)
output, error = p.communicate(input = test_string)

print output
print error
