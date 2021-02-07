import pytest
import numpy as np
import molellipsize as mes


def test_get_min_vol_ellipse(points):
    ET = mes.EllipsoidTool()

    test_center, test_radii, test_rotation = ET.get_min_vol_ellipse(
        points.points
    )

    print(points.center, points.radii, points.rotation)
    print(test_center, test_radii, test_rotation)

    for known, test_cent_dim in zip(points.center, test_center):
        assert np.isclose(known, test_cent_dim, rtol=0, atol=1E-1)

    for known, test_radius in zip(points.radii, test_radii):
        assert np.isclose(known, test_radius, rtol=0, atol=1E-2)

    assert np.allclose(
        points.rotation, test_rotation, rtol=0, atol=1E-3
    )
