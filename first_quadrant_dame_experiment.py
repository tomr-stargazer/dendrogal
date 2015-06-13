"""
This is an experiment to see if I can identify the structures associated 
with the Dame clouds DIRECTLY in the data.
Like, one-by-one.

THEN, identify their 'structure' properties and use that to re-derive the criteria for cloud definitions.
(Or at least justify it. I am pretty weak on that point these days.)

This script now also generates the paper figures involved & saves them to pdf.

"""

from __future__ import division
import os.path

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

from dendrogal.production.make_firstquad_stub import d
from dendrogal.production.cloud_extractor_q1 import first_quad_cloud_catalog

from dendrogal.integrated_viewer import IntegratedViewer
from dendrogal.production.integrated_map_figure import IntegratedMapFigure
from dendrogal.comparison_to_other_catalogs import plot_dame_ellipses_on_integrated_viewer, plot_dame_ellipses_on_imf
from dendrogal.production.dame_color_dict import dame_cmap
from dendrogal.production.integrated_map_figure import colorbrewer_blue

from dendrogal.hub_selection_manager import load_template
from astrodendro.scatter import Scatter

output_path = os.path.expanduser("~/Dropbox/Grad School/Research/Milkyway/paper/")


dv = d.viewer()

imf = IntegratedMapFigure(d.data, d.wcs, integration_limits=(-1.5, 1.5))

imf.ax.set_xlim(135, 496)
imf.ax.set_ylim(150, 336)

plot_dame_ellipses_on_imf(imf)

lon = imf.ax.coords['glon']
lon.set_ticks(spacing=5*u.deg, color='white', exclude_overlapping=True)
lon.display_minor_ticks(True)
lon.set_axislabel(r"$l$ (deg)", minpad=1.5)

vlsr = imf.ax.coords['vopt']
vlsr.set_ticks(spacing=25*u.m/u.s, color='white', exclude_overlapping=True) # erroneous units - why!?
vlsr.display_minor_ticks(True)
vlsr.set_axislabel(r"$v_{LSR}$ (km s$^{-1}$)")
vlsr.set_ticklabel_position('lr')

cube, header, selection_dictionary, selection_ID_dictionary = load_template('Dame86_clouds_selections', dv=dv, selection_key=2)

for [struct] in selection_dictionary.values():
    # the following code is cloned from the thumbnail thingy
    cloud_mask = d[struct.idx].get_mask()
    mask_lv = np.sum(cloud_mask, axis=1).astype('bool')

    imf.ax.contour(mask_lv, levels=[0.5], colors=colorbrewer_blue, linewidths=2, zorder=0.85)
    imf.ax.contour(mask_lv, levels=[0.5], colors='white', linewidths=3, zorder=0.8)

imf.fig.canvas.draw()
plt.show()
imf.fig.savefig(output_path+'Dame86_contours.pdf', bbox_inches='tight')

catalog = first_quad_cloud_catalog()
idx_list = [struct.idx for [struct] in selection_dictionary.values()]
subcat = catalog[ np.in1d(catalog['_idx'], idx_list) ]

plot_fig = plt.figure(figsize=(4,3))
ax = plot_fig.add_subplot(111)
ax.plot(subcat['v_rms'], subcat['n_descendants'], 'o', color=colorbrewer_blue)
ax.set_xlim(0,8)
ax.set_ylim(-0.2, 10.2)

ax.set_xlabel("$\sigma_v$ (km/s)")
ax.set_ylabel("n_descendants")

plot_fig.savefig(output_path+'Dame_clouds_vrms_ndesc.pdf', bbox_inches='tight')


