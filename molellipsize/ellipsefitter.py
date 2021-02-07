"""
EllisoidTool
============

From: https://github.com/minillinim/ellipsoid

ellipsoid.py
Copyright (C) 2018  Michael Imelfort

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

#. :class:`.EllipsoidTool`

EllipsoidTool class.

"""


from __future__ import division
import numpy as np


class EllipsoidTool:
    """
    EllipsoidTool for fitting ellipsoids.

    """

    def __init__(self):
        pass

    def get_min_vol_ellipse(
        self,
        points,
        tolerance=0.01,
    ):
        """
        Find the minimum volume ellipsoid which holds all the points.

        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/
        scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!

        Parameters
        ----------
        points : :class:`np.ndarray`
            (N, 3) numpy array of N points.

        tolerance : :class:`float`, optional
            Tolerance used in Khachiyan algorithm.

        Returns
        ------
        :class:`tuple` of :class:`np.ndarray`
            Tuple of arrays containing coordinates of ellipsoid center,
            ellipsoid radii, ellipsoid rotation matrix.

        """

        (N, d) = np.shape(points)
        d = float(d)

        # Q will be our working array
        Q = np.vstack([np.copy(points.T), np.ones(N)])
        QT = Q.T

        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        count = 0
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            # M the diagonal vector of an NxN matrix
            M = np.diag(np.dot(QT, np.dot(np.linalg.inv(V), Q)))
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0)
            step_size = step_size / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u
            count += 1

        # center of the ellipse
        center = np.dot(points.T, u)

        # the A matrix for the ellipse
        A = np.linalg.inv(
            np.dot(points.T, np.dot(np.diag(u), points)) -
            np.array([[a * b for b in center] for a in center])
        ) / d

        # Get the values we'd like to return
        U, s, rotation = np.linalg.svd(A)
        radii = 1.0/np.sqrt(s)

        return (center, radii, rotation)
