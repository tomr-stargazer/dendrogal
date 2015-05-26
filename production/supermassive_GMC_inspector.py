"""
This module creates thumbnails for the 10^7 Msun GMCs.

They're in the inner Galaxy (quadrants 1 and 2) and I bet some of them
have mis-assigned distance resolution.

"""

from __future__ import division

import os.path

import numpy as np
import matplotlib.pyplot as plt

from dendrogal.production.cloud_catalog_combiner import (extract_and_combine_catalogs, 
    export_firstquad_catalog, export_fourthquad_catalog, export_secondquad_catalog, 
    export_thirdquad_catalog, export_carina_catalog)
from dendrogal.production.cloud_extractor_q1 import d as quad1_d
from dendrogal.production.cloud_extractor_q2 import d as quad2_d
from dendrogal.production.cloud_extractor_q3 import d as quad3_d
from dendrogal.production.cloud_extractor_carina import d as carina_d
from dendrogal.production.cloud_extractor_q4 import d as quad4_d
from dendrogal.production.map_dendrogram_thumbnail_figure import make_thumbnail_dendro_figure

super_directory = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/supermassive/")
massive_directory = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/massive/")

allcat = extract_and_combine_catalogs()
cat1 = export_firstquad_catalog()
cat2 = export_secondquad_catalog()
cat3 = export_thirdquad_catalog()
cat_carina = export_carina_catalog()
cat4 = export_fourthquad_catalog()

supermassive_cat = allcat[allcat['mass'] > 5e6]
massive_cat = allcat[(allcat['mass'] > 1e6) & (allcat['mass'] <= 5e6)]

# print supermassive_cat['_idx']
# print supermassive_cat['survey']

survey_d_map = {8. : quad1_d, 
                17. : quad2_d,
                31. : quad3_d,
                33. : carina_d,
                36. : quad4_d}

survey_cat_map = {8. : cat1, 
                  17. : cat2,
                  31. : cat3,
                  33. : cat_carina,
                  36. : cat4}


for row in supermassive_cat:

    idx = row['_idx']
    d = survey_d_map[row['survey']]
    cat = survey_cat_map[row['survey']]

    fig = make_thumbnail_dendro_figure(d, cat, idx)
    catalog_properties = """
    D: {0} +{1} -{2} kpc
    KDA: {5}
    M: {3}x10^{4} Msun 
    """.format(row['distance'], row['error_distance_plus'], row['error_distance_minus'], str(row['mass'])[0], 
        int(np.floor(np.log10(row['mass']))), row['KDA_resolution'] )
    fig.text(0.2, -0.15, catalog_properties)

    fig.savefig(super_directory+"{0}.pdf".format(idx), bbox_inches='tight')
    print idx
    

for row in massive_cat:
    

    idx = row['_idx']
    d = survey_d_map[row['survey']]
    cat = survey_cat_map[row['survey']]

    fig = make_thumbnail_dendro_figure(d, cat, idx)

    catalog_properties = """
    D: {0} +{1} -{2} kpc
    KDA: {5}
    M: {3}x10^{4} Msun 
    """.format(row['distance'], row['error_distance_plus'], row['error_distance_minus'], str(row['mass'])[0], 
        int(np.floor(np.log10(row['mass']))), row['KDA_resolution'] )
    fig.text(0.2, -0.15, catalog_properties)
    

    fig.savefig(massive_directory+"{0}.pdf".format(idx), bbox_inches='tight')    

