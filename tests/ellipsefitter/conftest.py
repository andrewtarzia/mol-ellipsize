import pytest
import numpy as np
import molellipsize as mes

class CaseData:
    """
    A test case.

    Attributes
    ----------
    points : :class:`np.ndarray`
        (N, 3) array of points to fit ellipsoid to.

    center : :class:`list`
        Center of the ellipsoid.

    radii : :class:`tuple`
        Ellipsoid radii.

    rotation : :class:`np.ndarray`
        (3, 3) array of ellipsoid rotation matrix.

    """

    def __init__(self, points, center, radii, rotation):
        self.points = points
        self.center = center
        self.radii = radii
        self.rotation = rotation


@pytest.fixture(
    params=(
        CaseData(
            points=np.array([
                [1, 1, 1], [1, -1, 1], [-1, 1, 1], [-1, -1, 1],
                [1, 1, -1], [1, -1, -1], [-1, 1, -1], [-1, -1, -1],
            ]),
            center=[0, 0, 0],
            radii=(1.73205081, 1.73205081, 1.73205081),
            rotation=np.array([
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
            ]),
        ),
        CaseData(
            points=np.array([
                [1, 0, 0], [-1, 0, 0],
                [0, 1, 0], [0, -1, 0],
                [0, 0, 1], [0, 0, -1],
            ]),
            center=[0, 0, 0],
            radii=(1.0, 1.0, 1.0),
            rotation=np.array([
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
            ]),
        ),
        CaseData(
            points=np.array([
                [1, 0, 0], [-1, 0, 0],
                [0, 2, 0], [0, -2, 0],
                [0, 0, 1], [0, 0, -1],
            ]),
            center=[0, 0, 0],
            radii=(1.0, 1.0, 2.0),
            rotation=np.array([
                [0., 0., 1.],
                [1., 0., 0.],
                [0., 1., 0.],
            ]),
        ),
        CaseData(
            points=np.array([
                [2, 0, 0], [-2, 0, 0],
                [0, 1, 0], [0, -1, 0],
                [0, 0, 5], [0, 0, -5],
            ]),
            center=[0, 0, 0],
            radii=(1.0, 2.0, 5.0),
            rotation=np.array([
                [0., 1., 0.],
                [1., 0., 0.],
                [0., 0., 1.],
            ]),
        ),
        CaseData(
            points=np.array([
                [3, 1, 1], [1, 1, 1],
                [1, 2, 1], [1, 0, 1],
                [1, 1, 4], [1, 1, -2],
            ]),
            center=[1.5, 1.0, 1.0],
            radii=(1.05, 1.49, 3.13),
            rotation=np.array([
                [9.76950060e-03, -9.99952273e-01, -9.24355045e-05],
                [ 9.99943311e-01,  9.76980441e-03, -4.23378642e-03],
                [ 4.23448744e-03, -5.10682855e-05,  9.99991033e-01],
            ]),
        ),
    ),
)
def points(request):
    return request.param
