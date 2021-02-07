import pytest
from rdkit.Chem import AllChem as Chem
import molellipsize as mes


class CaseData:
    """
    A test case.

    Attributes
    ----------
    atom_atomic_nums : :class:`list`
        List of atom atomic numbers.

    pos_mat : :class:`iter` of :class:`list`
        List of atom coordinates.

    ellipse_center : :class:`tuple`
        Center of test molecule's ellipsoid.

    ellipse_diameters : :class:`tuple`
        Diameters of test molecule's ellipsoid.

    inertial_ratios : :class:`tuple`
        Inertial ratios of the test molecule.

    """

    def __init__(
        self,
        atom_atomic_nums,
        pos_mat,
        ellipse_center,
        ellipse_diameters,
        inertial_ratios
    ):

        # Make rdkit molecule.
        mol = Chem.EditableMol(Chem.Mol())
        for an in atom_atomic_nums:
            rdkit_atom = Chem.Atom(an)
            mol.AddAtom(rdkit_atom)

        mol = mol.GetMol()
        rdkit_conf = Chem.Conformer(len(atom_atomic_nums))
        for atom_id, atom_coord in enumerate(pos_mat):
            rdkit_conf.SetAtomPosition(atom_id, atom_coord)
            mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)
        mol.AddConformer(rdkit_conf)

        self.molecule = mol
        self.conformers = [0]
        self.ellipse_center = ellipse_center
        self.ellipse_diameters = ellipse_diameters
        self.inertial_ratios = inertial_ratios


@pytest.fixture(
    params=(
        # Oblate, in xy dir
        CaseData(
            atom_atomic_nums=[1, 1, 1, 1],
            pos_mat=[
                [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
            ],
            ellipse_center=(0.50, 0.50, 0.0),
            ellipse_diameters=(2.50904, 3.56057, 3.58628),
            inertial_ratios=(0.5, 0.5),
        ),
        # Oblate, but checks that the directionality plays no role,
        # by being in yz dir.
        CaseData(
            atom_atomic_nums=[1, 1, 1, 1],
            pos_mat=[
                [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
            ],
            ellipse_center=(0.0, 0.50, 0.50),
            ellipse_diameters=(2.50904, 3.56057, 3.58628),
            inertial_ratios=(0.5, 0.5),
        ),
        # Sphere.
        CaseData(
            atom_atomic_nums=[1, 1, 1, 1, 1, 1, 1, 1],
            pos_mat=[
                [1, 1, 1], [1, -1, 1], [-1, 1, 1], [-1, -1, 1],
                [1, 1, -1], [1, -1, -1], [-1, 1, -1], [-1, -1, -1],
            ],
            ellipse_center=(0.0, 0.0, 0.0),
            ellipse_diameters=(5.29935, 5.33006, 5.35154),
            inertial_ratios=(1.0, 1.0),
        ),
        # Prolate.
        CaseData(
            atom_atomic_nums=[1, 1, 1, 1],
            pos_mat=[
                [5, 1, 0], [5, -1, 0], [-5, 1, 0], [-5, -1, 0],
            ],
            ellipse_center=(0.0, 0.0, 0.0),
            ellipse_diameters=(3.35897, 5.18066, 17.32586),
            inertial_ratios=(0.0, 1.0),
        ),
    ),
)
def molecule(request):
    return request.param
