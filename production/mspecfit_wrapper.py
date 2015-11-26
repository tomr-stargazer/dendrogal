"""
Takes a given catalog and returns the mass spectrum fits.

"""

from __future__ import division
import os
import subprocess

import numpy as np

from dendrogal.production.config import idl_code_path, idl_executable


def prepare_stuff_for_mspecfit(catalog, name):
    """
    Makes an output suitable for IDL mspecfit stuff.

    """

    mass = catalog['mass']
    mass_err = np.sqrt(catalog['error_mass_plus']*catalog['error_mass_minus'])

    fname_mass = os.path.join(idl_code_path, name+'mass.txt')
    fname_err = os.path.join(idl_code_path, name+'err.txt')

    np.savetxt(fname_mass, mass)
    np.savetxt(fname_err, mass_err)


def send_mspec_command_to_idl(name, notrunc=1):

    p1 = subprocess.Popen(["echo",
        'read_files_noise_experiment(\'{0}\', {1})'.format(name, notrunc)],
        stdout=subprocess.PIPE)
    p2 = subprocess.Popen([idl_executable], 
        stdin=p1.stdout, cwd=idl_code_path, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = p2.communicate()[0]

    return output


def parse_idl_output(output):

    # ok so output is a string.
# ; OUTPUTS:
# ;   fit -- the parameters of the fit: [N_0, M_0, gamma] for a
# ;          truncated power law, [0.0 , N_0, gamma] for a power law.

    split_output = output.split("\n")

    N_0 = float(split_output[-5])
    M_0 = float(split_output[-4])
    gamma = float(split_output[-3])

    output_dict = {}
    output_dict['N_0'] = N_0
    output_dict['M_0'] = M_0
    output_dict['gamma'] = gamma

    return output_dict


def get_mspec_fit(catalog, name, notrunc=1):
    """ Combines the above helper functions. """

    # first - make the files
    prepare_stuff_for_mspecfit(catalog, name)

    # next - tell IDL to read them and spew out stdout noise
    output = send_mspec_command_to_idl(name, notrunc)

    # finally - parse it and return the values as a dict
    return parse_idl_output(output)

