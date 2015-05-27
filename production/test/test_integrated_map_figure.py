"""
To be run with py.test.

"""

from __future__ import division

import numpy as np
from numpy.testing import assert_equal

from ..integrated_map_figure import sanitize_integration_limits


def test_sanitize_integration_limits():

    test_datacube = np.zeros((4,5,6))

    test_limits_ax0 = (-1, 3)
    test_limits_ax1 = (3, 8)
    test_limits_ax2 = (2, 5)

    san_limits_ax0 = sanitize_integration_limits(test_limits_ax0, test_datacube, axis=0)
    san_limits_ax1 = sanitize_integration_limits(test_limits_ax1, test_datacube, axis=1)
    san_limits_ax2 = sanitize_integration_limits(test_limits_ax2, test_datacube, axis=2)

    expected_ax0 = (0, 3)
    expected_ax1 = (3, 5)
    expected_ax2 = (2, 5)

    assert_equal(san_limits_ax0, expected_ax0)
    assert_equal(san_limits_ax1, expected_ax1)
    assert_equal(san_limits_ax2, expected_ax2)




