"""
Takes a given catalog and returns the mass spectrum fits.

"""

def multipanel_write_masses_and_errors_to_file_output():
    """
    Makes an output suitable for IDL mspecfit stuff.

    """

    sky_text_dict = {'north': 'Northern', 'south': 'Southern', 'north_plus_south' : 'All'}
    gal_text_dict = {'inner': ' inner Galaxy', 'outer': ' outer Galaxy', 'inner_plus_outer' : ' Galaxy (combined)'}

    sky_subset_dict = {}
    sky_subset_dict['north'] = pruned_catalog['x_cen'] < 180
    sky_subset_dict['south'] = pruned_catalog['x_cen'] > 180
    sky_subset_dict['north_plus_south'] = pruned_catalog['x_cen'] != 0

    gal_subset_dict = {}
    gal_subset_dict['inner'] = pruned_catalog['R_gal'] < 8
    gal_subset_dict['outer'] = pruned_catalog['R_gal'] > 9
    gal_subset_dict['inner_plus_outer'] = pruned_catalog['R_gal'] != 0

    # the core loop.
    for i, sky_region in enumerate(['north', 'south', 'north_plus_south']):

        for j, galactic_region in enumerate(['inner', 'outer', 'inner_plus_outer']):

            if (galactic_region == 'inner_plus_outer') and (sky_region != 'north_plus_south'):
                continue

            if galactic_region == 'outer':
                max_mass = 3e7
                min_mass = 1e4
            else:
                max_mass = 3e7
                min_mass = 1e5

            catalog = pruned_catalog[sky_subset_dict[sky_region] & gal_subset_dict[galactic_region] & 
                                    (~np.isnan(pruned_catalog['mass'])) &
                                    (pruned_catalog['mass']>min_mass) & (pruned_catalog['mass']<max_mass)]

            mass = catalog['mass']
            mass_err = np.sqrt(catalog['error_mass_plus']*catalog['error_mass_minus'])

            path = os.path.expanduser("~/Documents/Code/idl-low-sky/eroslib/")

            fname_mass = path + galactic_region+sky_region+'mass.txt'
            fname_err = path + galactic_region+sky_region+'err.txt'

            np.savetxt(fname_mass, mass)
            np.savetxt(fname_err, mass_err)

            print len(mass), sky_region, galactic_region