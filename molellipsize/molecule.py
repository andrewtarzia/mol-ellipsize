"""
Molecule
========

#. :class:`.Molecule`

Molecule class.

"""

from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Descriptors3D import NPR1, NPR2
from rdkit.Geometry import rdGeometry

import numpy as np

from .ellipsefitter import EllipsoidTool
from .utilities import plot_ellipsoid


class Molecule:
    """
    Molecule to calculate size of.

    """

    def __init__(self, rdkitmol, conformers):
        """
        Initialize a :class:`Molecule` instance.

        Parameters
        ----------
        rdkitmol : :class:`RDKit.Molecule`
            RDKit molecule to get size of.

        conformers : :class:`iterable`
            Iterable of the conformer ids used to access the conformers
            in the rdkit molecule.

        """

        self._rdkitmol = rdkitmol
        self._conformers = conformers

    def get_inertial_ratios(self):
        """
        Get inertial 3D descriptors for all conformers in mol.

        Returns
        -------
        conf_ratios : :class:`dict` of :class:`tuple`
            Dictionary of ratio_1 and ratio_2 of all conformers.
            Key is conformer id. Ratio 1 is I1/I3, ratio 2 is I2/I3.

        """

        conf_ratios = {}
        for cid in self._conformers:
            conf_ratios[cid] = (
                NPR1(self._rdkitmol, confId=cid),
                NPR2(self._rdkitmol, confId=cid),
            )

        return conf_ratios

    def get_molecule_shape(
        self,
        conformer,
        cid,
        vdwscale,
        boxmargin,
        spacing
    ):
        """
        Get the shape of a conformer of a molecule as a grid.

        """
        box = Chem.ComputeConfBox(conformer)
        sideLen = (
            box[1].x-box[0].x + 2*boxmargin,
            box[1].y-box[0].y + 2*boxmargin,
            box[1].z-box[0].z + 2*boxmargin,
        )
        shape = rdGeometry.UniformGrid3D(
            2*sideLen[0],
            2*sideLen[1],
            2*sideLen[2],
            spacing=spacing
        )
        Chem.EncodeShape(
            self._rdkitmol,
            shape,
            confId=cid,
            ignoreHs=False,
            vdwScale=vdwscale
        )
        return box, sideLen, shape

    def get_hitpoints(self, shape):
        """
        Get points with value > 2 that are within vdw shape.

        """

        hit_points = []
        for idx in range(shape.GetSize()):
            value = shape.GetVal(idx)
            if value > 2:
                pt = shape.GetGridPointLoc(idx)
                point = np.array([pt.x, pt.y, pt.z])
                hit_points.append(point)
        hit_points = np.asarray(hit_points)

        return hit_points

    def get_ellipsoids(
        self,
        vdwscale,
        boxmargin,
        spacing,
    ):
        """
        Get min volume ellipsoids for all conformers in mol.

        Good values:
            vdwScale:0.9
            boxMargin:4.0
            spacing:0.5

        Returns
        -------
        conf_ellipsoids : :class:`dict` of :class:`tuple`
            Dictionary of conformer ellipsoids, key is conformer id.
            Value is (center, diameters, rotation matrix).

        """

        conf_ellipsoids = {}
        for cid in self._conformers:
            conformer = self._rdkitmol.GetConformer(cid)
            box, sideLen, shape = self.get_molecule_shape(
                conformer=conformer,
                cid=cid,
                vdwscale=vdwscale,
                boxmargin=boxmargin,
                spacing=spacing,
            )

            hit_points = self.get_hitpoints(shape)

            # Find the ellipsoid that envelopes all hit points.
            ET = EllipsoidTool()
            (center, radii, rotation) = ET.get_min_vol_ellipse(
                points=hit_points,
                tolerance=0.01,
            )
            diameters = list(np.sort(radii)*2)

            conf_ellipsoids[cid] = (center, diameters, rotation)

        return conf_ellipsoids

    def draw_conformer_hitpoints(
        self,
        cid,
        center,
        diameter,
        rotation,
        vdwscale,
        boxmargin,
        spacing,
        filename,
    ):
        """
        Draw hitpoints of ellipsoid fit in plot.

        """

        conformer = self._rdkitmol.GetConformer(cid)
        box, sideLen, shape = self.get_molecule_shape(
            conformer=conformer,
            cid=cid,
            vdwscale=vdwscale,
            boxmargin=boxmargin,
            spacing=spacing,
        )

        hit_points = self.get_hitpoints(shape)
        fig, ax = plot_ellipsoid(
            center,
            diameter,
            rotation,
        )
        ax.scatter(
            hit_points[:, 0], hit_points[:, 1], hit_points[:, 2],
            color='g',
            marker='x',
            edgecolor=None,
            s=50,
            alpha=0.5,
        )

        fig.tight_layout()
        fig.savefig(filename, dpi=720, bbox_inches='tight')

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'<{self.__class__.__name__} at {id(self)}> '
            f'with {len(self._conformers)} conformers'
        )
