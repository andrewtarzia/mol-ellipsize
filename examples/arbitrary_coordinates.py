from rdkit.Chem import AllChem as Chem
import numpy as np
import molellipsize as mes


def main():

    # Using a Molecule.
    # Fake the positions with C atoms, and some fake bonding.
    atom_atomic_nums = [6, 6, 6, 6, 6, 6]
    # Atom1 id, atom2 id, bond order.
    bonds = [(0, 1, 1), (3, 4, 1)]
    atom_positions = [
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1],
    ]
    # Make rdkit molecule.
    mol = Chem.EditableMol(Chem.Mol())
    for an in atom_atomic_nums:
        rdkit_atom = Chem.Atom(an)
        mol.AddAtom(rdkit_atom)

    for bond in bonds:
        mol.AddBond(
            beginAtomIdx=bond[0],
            endAtomIdx=bond[1],
            order=Chem.BondType(bond[2]),
        )

    mol = mol.GetMol()
    rdkit_conf = Chem.Conformer(len(atom_atomic_nums))
    for atom_id, atom_coord in enumerate(atom_positions):
        rdkit_conf.SetAtomPosition(atom_id, atom_coord)
        mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)
    mol.AddConformer(rdkit_conf)
    # Hard set this.
    conformers = [0]
    mes_mol = mes.Molecule(mol, conformers)
    conf_ratios = mes_mol.get_inertial_ratios()
    conf_ellipsoids = mes_mol.get_ellipsoids(1.0, 4.0, 0.5)
    center = conf_ellipsoids[0][0]
    diameters = conf_ellipsoids[0][1]
    print(
        f'This fake molecule makes a nice sphere at {center}'
        f' with diameters: {diameters} and inertial ratios '
        f'{conf_ratios[0]}'
    )

    # Without a molecule.
    points = np.array([
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1],
    ])
    ET = mes.EllipsoidTool()
    center, radii, rotation = ET.get_min_vol_ellipse(points)
    print(
        f'These points make a nice sphere at {center}'
        f' with radii: {radii}'
    )


if __name__ == '__main__':
    main()
