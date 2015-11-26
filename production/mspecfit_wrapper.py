"""
Takes a given catalog and returns the mass spectrum fits.

"""

from __future__ import division
import os
import subprocess

import numpy as np

from dendrogal.production.config import idl_code_path


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

    cwd = os.getcwd()

    idl_command = "/Applications/exelis/idl83/bin/idl"


    # os.chdir(idl_code_path)

    p1 = subprocess.Popen(["echo", 'read_files_noise_experiment(\'{0}\', {1})'.format(name, notrunc)], stdout=subprocess.PIPE)

    # return p1.communicate()

    p2 = subprocess.Popen([idl_command], stdin=p1.stdout, cwd=idl_code_path, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # p2 = subprocess.Popen([idl_command], stdin='\"read_files_noise_experiment(\'{0}\', {1})\"'.format(name, notrunc), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # p1.communicate()
    output = p2.communicate()[0]

    # os.chdir(cwd)


    # output = os.system("echo \"read_files_noise_experiment(\'{0}\', {1})\" | {2}".format(name, notrunc, idl_command))

    return output


def send_mspec_command_to_idl_backup(name, notrunc=1):

    idl_command = "/Applications/exelis/idl83/bin/idl"

    cwd = os.getcwd()

    # p1 = subprocess.Popen(["echo", '\"read_files_noise_experiment(\'{0}\', {1})\"'.format(name, notrunc)], stdout=subprocess.PIPE)

    # # return p1.communicate()


    # p2 = subprocess.Popen([idl_command], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=idl_code_path)

    # p2 = subprocess.Popen([idl_command], stdin='\"read_files_noise_experiment(\'{0}\', {1})\"'.format(name, notrunc), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # p1.communicate()
    # output = p2.communicate()

    os.chdir(idl_code_path)

    output = os.system("echo \"read_files_noise_experiment(\'{0}\', {1})\" | {2}".format(name, notrunc, idl_command))

    os.chdir(cwd)

    return output


def send_mspec_command_to_idl_backup2(name, notrunc=1):

    idl_command = "/Applications/exelis/idl83/bin/idl"

    p1 = subprocess.check_output("echo \"read_files_noise_experiment(\'{0}\', {1})\" | {2}".format(name, notrunc, idl_command))

    return p1

    # p1 = subprocess.Popen(["echo", '\"read_files_noise_experiment(\'{0}\', {1})\"'.format(name, notrunc)], stdout=subprocess.PIPE)

    # # return p1.communicate()


    # p2 = subprocess.Popen([idl_command], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=idl_code_path)

    # p2 = subprocess.Popen([idl_command], stdin='\"read_files_noise_experiment(\'{0}\', {1})\"'.format(name, notrunc), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # p1.communicate()
    # output = p2.communicate()

    os.chdir(idl_code_path)

    output = os.system("echo \"read_files_noise_experiment(\'{0}\', {1})\" | {2}".format(name, notrunc, idl_command))

    os.chdir(cwd)

    return output    


def parse_idl_output():
    pass
