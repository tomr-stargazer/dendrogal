"""
A custom colormap based on RgGy_r that is "stretched" in a way I think Tom Dame will like.

"""

from __future__ import division

import numpy as np

# custom colormap - inspired by
# http://stackoverflow.com/questions/16152052/matplotlib-python-change-single-color-in-colormap/16163481
from matplotlib.colors import LinearSegmentedColormap

halfway = 0.6
first_half_power = 1/3
second_half_power = 3/4

stretches_first_half = (np.linspace(0, halfway, 6) / halfway)**(first_half_power) * halfway
xx, stretch_1, stretch_2, stretch_3, stretch_4, stretch_5 = stretches_first_half

stretches_second_half = (np.linspace(halfway, 1, 6) / halfway)**(second_half_power) * halfway
xx, stretch_6, stretch_7, stretch_8, stretch_9, yy = stretches_second_half

# stretch_1 = 0.1
# stretch_2 = 0.2
# stretch_3 = 0.3
# stretch_4 = 0.4
# stretch_5 = 0.5
# stretch_6 = 0.6
# stretch_7 = 0.7
# stretch_8 = 0.8
# stretch_9 = 0.9

# derived from:
# from matplotlib import cm
# cm.RdGy_r._segmentdata

color_dict = \
{u'blue': [(0.0, 0.10196078568696976, 0.10196078568696976),
  (stretch_1, 0.3019607961177826, 0.3019607961177826),
  (stretch_2, 0.529411792755127, 0.529411792755127),
  (stretch_3, 0.729411780834198, 0.729411780834198),
  (stretch_4, 0.8784313797950745, 0.8784313797950745),
  (stretch_5, 1.0, 1.0),
  (stretch_6, 0.7803921699523926, 0.7803921699523926),
  (stretch_7, 0.5098039507865906, 0.5098039507865906),
  (stretch_8, 0.3019607961177826, 0.3019607961177826),
  (stretch_9, 0.16862745583057404, 0.16862745583057404),
  (1.0, 0.12156862765550613, 0.12156862765550613)],
 u'green': [(0.0, 0.10196078568696976, 0.10196078568696976),
  (stretch_1, 0.3019607961177826, 0.3019607961177826),
  (stretch_2, 0.529411792755127, 0.529411792755127),
  (stretch_3, 0.729411780834198, 0.729411780834198),
  (stretch_4, 0.8784313797950745, 0.8784313797950745),
  (stretch_5, 1.0, 1.0),
  (stretch_6, 0.8588235378265381, 0.8588235378265381),
  (stretch_7, 0.6470588445663452, 0.6470588445663452),
  (stretch_8, 0.3764705955982208, 0.3764705955982208),
  (stretch_9, 0.0941176488995552, 0.0941176488995552),
  (1.0, 0.0, 0.0)],
 u'red': [(0.0, 0.10196078568696976, 0.10196078568696976),
  (stretch_1, 0.3019607961177826, 0.3019607961177826),
  (stretch_2, 0.529411792755127, 0.529411792755127),
  (stretch_3, 0.729411780834198, 0.729411780834198),
  (stretch_4, 0.8784313797950745, 0.8784313797950745),
  (stretch_5, 1.0, 1.0),
  (stretch_6, 0.9921568632125854, 0.9921568632125854),
  (stretch_7, 0.95686274766922, 0.95686274766922),
  (stretch_8, 0.8392156958580017, 0.8392156958580017),
  (stretch_9, 0.6980392336845398, 0.6980392336845398),
  (1.0, 0.40392157435417175, 0.40392157435417175)]}

dame_cmap = LinearSegmentedColormap('dame_cmap', color_dict)
