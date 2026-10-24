"""Terrain controls (section 7.4) on synthetic lines and rasters."""

import math

import numpy as np
import pytest
from rasterio.transform import from_origin
from shapely.geometry import LineString

from floodforest.controls import curvature_per_km, gradient_per_km, max_within, outer_bend


def test_gradient_is_the_local_fall_per_km():
    km = np.arange(0, 3.01, 0.1)
    z0 = 2000 - 20 * km  # 20 m per km
    z0[5] = np.nan
    grad = gradient_per_km(km, z0, 500)
    assert np.allclose(grad, -20)


def test_curvature_sign_and_outer_bank():
    # East for 1 km, then a left turn to the north.
    river = LineString([(0, 0), (1000, 0), (1000, 1000)])
    curv = curvature_per_km(river, np.array([0.0, 500.0, 1000.0]), 300)
    assert math.isnan(curv[0])
    assert curv[1] == pytest.approx(0.0)
    assert curv[2] == pytest.approx((math.pi / 2) / 0.6)
    bend = outer_bend(np.array([curv[2], curv[2], -1.0, 0.0]),
                      np.array(["left", "right", "left", "right"]))
    assert bend.tolist() == [0, 1, 1, 0]


def test_max_within_uses_a_circle():
    values = np.zeros((101, 101))
    values[50, 80] = 5.0   # 300 m east of the centre pixel
    values[80, 80] = 9.0   # about 424 m away, outside 300 m
    transform = from_origin(0, 1010, 10, 10)
    centre = np.array([[505.0, 505.0]])
    assert max_within(values, transform, centre, 300)[0] == 5.0
