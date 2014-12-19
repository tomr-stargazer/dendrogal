def disqualify_and_return_copy(input_catalog):

    catalog = input_catalog.copy(copy_data=True)

    disqualified = ((catalog['mass'] > 1e7) | 
                    (catalog['mass'] < 5e3) | 
                    (catalog['major_sigma'] > 10) | 
                    (catalog['v_rms'] > 30) | 
                    (catalog['size'] > 1000) |
                    # (np.abs(catalog['v_cen']) < 13) |
                    # (np.abs(catalog['x_cen'] - 180) < 10) |
                    (catalog['area_exact'] > 50) | 
                    (catalog['radius'] < 0.18) | 
                    (catalog['radius'] > 2*1.32))

    catalog['Distance'][disqualified] = np.nan
    catalog['v_rms'][disqualified] = np.nan
    catalog['size'][disqualified] = np.nan
    catalog['virial'][disqualified] = np.nan
    catalog['mass'][disqualified] = np.nan

    return catalog