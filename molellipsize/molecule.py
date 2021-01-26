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


class Molecule:
    """
    Molecule to calculate size of.

    """

    def __init__(self, rdkitmol, conformers):
        """
        Initialize a :class:`Bond` instance.

        Parameters
        ----------
        rdkitmol : :class:`RDKit.Molecule`
            RDKit molecule to get size of.

        conformers : :class:`None`
            X

        """

        self._rdkitmol = rdkitmol
        self._conformers = conformers

    def get_inertial_ratios(self):
        """
        Get inertial 3D descriptors for all conformers in mol.

        Returns
        ------
        XXXX
        :class:`tuple` of :class:`float`
            Yields tuple of conformer id, ratio_1 and ratio_2.
            Ratio 1 is I1/I3, ratio 2 is I2/I3.

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

    def get_diameters(self, vdwscale, boxmargin, spacing):
        """
        Get XX for all conformers in mol.

        Good values:
            vdwScale:0.9
            boxMargin:4.0
            spacing:0.5

        Returns
        -------
        :class:`tuple` of :class:`float`
            Yields tuple of conformer id, [d1, d2, d3]

        """

        conf_diameters = {}
        for cid in self._conformers:
            conformer = self._rdkitmol.GetConformer(cid)
            box, sideLen, shape = self.get_molecule_shape(
                conformer=conformer,
                cid=cid,
                vdwscale=vdwscale,
                boxmargin=boxmargin,
                spacing=spacing,
            )

            # Get ellipsoid fitting all points with value > 2.
            # - i.e. within vdw shape.
            hit_points = []
            for idx in range(shape.GetSize()):
                pt = shape.GetGridPointLoc(idx)
                value = shape.GetVal(idx)
                if value > 2:
                    point = np.array([pt.x, pt.y, pt.z])
                    hit_points.append(point)
            hit_points = np.asarray(hit_points)

            # Find the ellipsoid that envelopes all hit points.
            ET = EllipsoidTool()
            (center, radii, rotation) = ET.getMinVolEllipse(
                P=hit_points,
                tolerance=0.01,
                do_step_plot=False,
            )

            conf_diameters[cid] = sorted(np.asarray(radii)*2)

        return conf_diameters

    def draw_hitpoints(self):
        """
        Draw hitpoints of ellipsoid fit in plot.

        """
        pass
        hit_point_plot(
            hit_points,
            ET,
            center,
            radii,
            rotation
        )
