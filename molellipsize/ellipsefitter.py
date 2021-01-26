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
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg


class EllipsoidTool:
    """
    EllipsoidTool for fitting ellipsoids.

    """

    def __init__(self):
        pass

    def step_plot(self, step, err, tolerance, u, d, P):
        print(step, err, tolerance)
        # center of the ellipse
        center = np.dot(P.T, u)

        # the A matrix for the ellipse
        A = linalg.inv(
            np.dot(P.T, np.dot(np.diag(u), P)) -
            np.array([[a * b for b in center] for a in center])
        ) / d

        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)
        print(radii)
        print(self.getEllipsoidVolume(radii=radii))

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # plot points
        ax.scatter(
            P[:, 0], P[:, 1], P[:, 2],
            color='#CB4335', marker='o', s=20,
            alpha=0.4
        )

        # plot ellipsoid
        self.plotEllipsoid(
            center, radii, rotation, ax=ax, plotAxes=True,
            cageAlpha=0.6
        )
        plt.axis('off')

        fig.tight_layout()
        fig.savefig(
            f"ellips_plots/ellip_step_{step}.pdf",
            dpi=720,
            bbox_inches='tight'
        )
        plt.close(fig)
        del fig
        input()

    def getMinVolEllipse(
        self,
        P=None,
        tolerance=0.01,
        do_step_plot=False
    ):
        """ Find the minimum volume ellipsoid which holds all the points

        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/
        scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!

        Here, P is a numpy array of N dimensional points like this:
        P = [[x,y,z,...], <-- one point per line
             [x,y,z,...],
             [x,y,z,...]]

        Returns:
        (center, radii, rotation)

        """
        (N, d) = np.shape(P)
        d = float(d)

        # Q will be our working array
        Q = np.vstack([np.copy(P.T), np.ones(N)])
        QT = Q.T

        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        count = 0
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            # M the diagonal vector of an NxN matrix
            M = np.diag(np.dot(QT, np.dot(linalg.inv(V), Q)))
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0)
            step_size = step_size / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u
            count += 1
            if do_step_plot:
                self.step_plot(count, err, tolerance, u, d, P)

        # center of the ellipse
        center = np.dot(P.T, u)

        # the A matrix for the ellipse
        A = linalg.inv(
            np.dot(P.T, np.dot(np.diag(u), P)) -
            np.array([[a * b for b in center] for a in center])
        ) / d

        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)

        return (center, radii, rotation)

    def getEllipsoidVolume(self, radii):
        """Calculate the volume of the blob"""
        return 4./3.*np.pi*radii[0]*radii[1]*radii[2]

    def plotEllipsoid(
        self,
        center,
        radii,
        rotation,
        ax=None,
        plotAxes=False,
        cageColor='k',
        cageAlpha=0.2
    ):
        """Plot an ellipsoid"""
        make_ax = ax is None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)

        # cartesian coordinates that correspond to the spherical
        # angles:
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i, j], y[i, j], z[i, j]] = np.dot(
                    [x[i, j], y[i, j], z[i, j]],
                    rotation
                ) + center

        if plotAxes:
            # make some purdy axes
            axes = np.array([[radii[0], 0.0, 0.0],
                             [0.0, radii[1], 0.0],
                             [0.0, 0.0, radii[2]]])
            # rotate accordingly
            for i in range(len(axes)):
                axes[i] = np.dot(axes[i], rotation)

            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                Z3 = np.linspace(-p[2], p[2], 100) + center[2]
                ax.plot(X3, Y3, Z3, color=cageColor)

        # plot ellipsoid
        ax.plot_wireframe(
            x, y, z,
            rstride=4,
            cstride=4,
            color=cageColor,
            alpha=cageAlpha
        )

        ax.set_xlabel(r'$x$ [$\mathrm{\AA}$]')
        ax.set_ylabel(r'$y$ [$\mathrm{\AA}$]')
        ax.set_zlabel(r'$z$ [$\mathrm{\AA}$]')

        # along
        # ax.view_init(0, 0)
        # across
        # ax.view_init(90, 0)

        if make_ax:
            plt.show()
            plt.close(fig)
            del fig
