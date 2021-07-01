"""
A script to compare calculed molecular size with rdkit 3D descriptors.

Uses the platinum database
(https://comp3d.univie.ac.at/servers-datasets/#c523270)

"""

import json
from os import mkdir
from os.path import exists
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolDescriptors
from scipy.stats import spearmanr
import pandas as pd
import matplotlib.pyplot as plt

import molellipsize as mes


def desc_list():
    return [
        rdMolDescriptors.CalcPBF, rdMolDescriptors.CalcPMI1,
        rdMolDescriptors.CalcPMI2, rdMolDescriptors.CalcPMI3,
        rdMolDescriptors.CalcNPR1, rdMolDescriptors.CalcNPR2,
        rdMolDescriptors.CalcRadiusOfGyration,
        rdMolDescriptors.CalcInertialShapeFactor,
        rdMolDescriptors.CalcAsphericity,
        rdMolDescriptors.CalcEccentricity,
        rdMolDescriptors.CalcSpherocityIndex,
    ]

def calculate_rdkit_descriptors(molecule, cid):
    molecule.ClearComputedProps()
    descriptor_functions = desc_list()

    return {
        i.__name__.replace('Calc', ''): i(molecule, confId=cid)
        for i in descriptor_functions
    }


def calculate_coeffs(data):

    coeffs = {}
    midellips = [
        data['diameters'][cid][1]
        for cid in data['diameters']
    ]
    rg_vs_size = []
    pmi_vs_size = []
    for desc in desc_list():
        name = desc.__name__.replace('Calc', '')
        values = [
            data['descriptors'][cid][name]
            for cid in data['descriptors']
        ]
        coeffs[name], _ = spearmanr(midellips, values)
        if name == 'RadiusOfGyration':
            for i, j in zip(midellips, values):
                rg_vs_size.append((i, j))
        if name == 'PMI2':
            for i, j in zip(midellips, values):
                pmi_vs_size.append((i, j))
    return coeffs, rg_vs_size, pmi_vs_size


def plot_all_rg_vs_molsize(data):

    print(data)
    fig, ax = plt.subplots(figsize=(5, 5))
    dx = [i[0] for i in data]
    dy = [i[1] for i in data]
    ax.scatter(
        dx, dy,
        c='gray',
        edgecolors='none',
        alpha=1.0,
        s=10
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'mol-size [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(r'rdkit $R_g$ [$\mathrm{\AA}$]', fontsize=16)

    fig.tight_layout()
    fig.savefig('rg_vs_mid.pdf', dpi=160, bbox_inches='tight')


def plot_all_pmi_vs_molsize(data):

    print(data)
    fig, ax = plt.subplots(figsize=(5, 5))
    dx = [i[0] for i in data]
    dy = [i[1] for i in data]
    ax.scatter(
        dx, dy,
        c='gray',
        edgecolors='none',
        alpha=1.0,
        s=10
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'mol-size [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(r'rdkit PMI2', fontsize=16)

    fig.tight_layout()
    fig.savefig('pmi_vs_mid.pdf', dpi=160, bbox_inches='tight')


def main():
    if not exists('rdkit_cf_data'):
        mkdir('rdkit_cf_data')

    # Read in Platnium database.
    molecules = [
        x for x in Chem.SDMolSupplier(
            fileName=(
                '/home/atarzia/Downloads/temp/platinum_comparison/'
                'platinum_dataset_2017_01.sdf'
            ),
            removeHs=False
        )
    ]
    # Assign atomic chirality based on the structures:
    for m in molecules:
        Chem.AssignAtomChiralTagsFromStructure(m)
    print(len(molecules))

    N = 50
    conformer_spearmen_coefficients = {}
    all_rg_molsize_pairs = []
    all_pmi_molsize_pairs = []
    for i, mol in enumerate(molecules[:700]):
        mol_data_file = f'rdkit_cf_data/mol_{i}.json'
        if exists(mol_data_file):
            with open(mol_data_file, 'r') as f:
                mol_data = json.load(f)
        else:
            mol_data = {}
            rdkitmol, conformers = mes.ETKDG_UFF_conformers(
                rdkitmol=mol,
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

            mol_data['diameters'] = diameters
            descriptors = {
                cid: calculate_rdkit_descriptors(rdkitmol, cid)
                for cid in conformers
            }
            mol_data['descriptors'] = descriptors
            with open(mol_data_file, 'w') as f:
                json.dump(mol_data, f)
        conformer_spearmen_coefficients[i], rg_vs_size, pmi_vs_size = (
            calculate_coeffs(mol_data)
        )
        for i in rg_vs_size:
            all_rg_molsize_pairs.append(i)
        for i in pmi_vs_size:
            all_pmi_molsize_pairs.append(i)

    df = pd.DataFrame.from_dict(conformer_spearmen_coefficients)
    df.to_csv('conf_spearmans.csv')
    plot_all_rg_vs_molsize(all_rg_molsize_pairs)
    plot_all_pmi_vs_molsize(all_pmi_molsize_pairs)


if __name__ == '__main__':
    main()
