from rdkit.Chem import AllChem as Chem
import molellipsize as mes


def main():

    dict_of_smiles = {
        'c6h6': 'c1ccccc1',
        'c6h12': 'C1CCCCC1',
        'c8': 'CCCCCCCC',
        'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    }

    N = 50

    for name in dict_of_smiles:
        smiles = dict_of_smiles[name]
        print(name, smiles)
        rdkitmol = Chem.MolFromSmiles(smiles)
        rdkitmol = Chem.AddHs(rdkitmol)
        Chem.SanitizeMol(rdkitmol)
        rdkitmol, conformers = mes.ETKDG_UFF_conformers(
            rdkitmol=rdkitmol,
            num_conformers=N,
            randomseed=1000,
        )
        mes_mol = mes.Molecule(
            rdkitmol=rdkitmol,
            conformers=conformers,
        )
        ellipsoids = mes_mol.get_ellipsoids(
            boxmargin=4.0,
            vdwscale=0.9,
            spacing=0.5,
        )
        diameters = {i: ellipsoids[i][1] for i in ellipsoids}
        ratios = mes_mol.get_inertial_ratios()
        mes.plot_shapes(ratios, filename=f'{name}_shape.pdf')
        mes.plot_diameters(diameters, filename=f'{name}_size.pdf')

        # Pick an ellipsoid to draw.
        cid = 0
        center, diameter, rotation = ellipsoids[cid]
        mes_mol.draw_conformer_hitpoints(
            cid=cid,
            center=center,
            diameter=diameter,
            rotation=rotation,
            boxmargin=4.0,
            vdwscale=0.9,
            spacing=0.5,
            filename=f'{name}_ellhp.pdf',
        )
        fig, ax = mes.plot_ellipsoid(
            center=center,
            diameter=diameter,
            rotation=rotation,
        )
        fig.tight_layout()
        fig.savefig(f'{name}_ell.pdf', dpi=720, bbox_inches='tight')


if __name__ == '__main__':
    main()
