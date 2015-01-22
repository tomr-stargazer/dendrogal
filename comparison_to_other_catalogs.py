"""
A script allowing (visual, statistical) comparisons between our catalog and others.

Some included here:
Dame+1986 (First Quadrant, 33 objects)
Solomon+1987 (First Quadrant, 464 objects)
Scoville+1987 (First Quadrant, 1427 objects)


"""

from __future__ import division

import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from astropy import table

data_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/other_catalogs/")

dame86_catalog = table.Table.read(data_path+'Dame+86_catalog_sanitized.txt', format='ascii', header_start=2)

id_column = dame86_catalog['ID']
vcen_column = [x.split(',')[1] for x in id_column]
lcen_column = dame86_catalog['LII']
sky_radius_column = np.degrees(dame86_catalog['R']/(dame86_catalog['D']*1000))

fig1 = plt.figure()

plt.plot(lcen_column, vcen_column, 'bo')
plt.gca().invert_xaxis()

plt.show()

fig2 = plt.figure()

ax = fig2.add_subplot(111)

ells = [Ellipse(xy=zip(lcen_column, vcen_column)[i], 
                width=2*sky_radius_column[i], height=dame86_catalog['DV'][i]) for i in range(len(dame86_catalog))]
for e in ells:
    ax.add_artist(e)
    e.set_facecolor('none')
plt.xlim(65,0)
plt.ylim(0,150)
plt.show()