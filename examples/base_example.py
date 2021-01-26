from rdkit.Chem import AllChem as Chem
import molellipsize as mes


def main():

    dict_of_smiles = {
        'c6h6': 'c1ccccc1',
        'c6h12': 'C1CCCCC1',
    }

    for name in dict_of_smiles:
        smiles = dict_of_smiles[name]
        print(name, smiles)
        rdkitmol = Chem.MolFromSmiles(smiles)
        rdkitmol = Chem.AddHs(rdkitmol)
        Chem.SanitizeMol(rdkitmol)
        rdkitmol, conformers = mes.ETKDG_UFF_conformers(
            rdkitmol=rdkitmol,
            num_conformers=10,
            randomseed=1000,
        )
        print(name, [i for i in conformers])
        mes_mol = mes.Molecule(
            rdkitmol=rdkitmol,
            conformers=conformers,
        )
        diameters = mes_mol.get_diameters(
            boxmargin=4.0,
            vdwscale=0.9,
            spacing=0.5,
        )
        ratios = mes_mol.get_inertial_ratios()
        mes.plot_shapes(ratios, filename=f'{name}_shape.pdf')
        mes.plot_diameters(diameters, filename=f'{name}_size.pdf')


if __name__ == '__main__':
    main()
