"""
This module defines general-purpose objects, functions and classes.

"""

import numpy as np
from rdkit.Chem import AllChem as Chem
import matplotlib.pyplot as plt


def ETKDG_UFF_conformers(rdkitmol, num_conformers, randomseed):
    """
    Generate conformers with RDKit.

    """
    # Use a set randomSeed so that running the code multiple times
    # gives the same series of conformers.
    conformers = Chem.EmbedMultipleConfs(
        mol=rdkitmol,
        numConfs=num_conformers,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True,
        randomSeed=randomseed
    )

    for cid in conformers:
        Chem.UFFOptimizeMolecule(rdkitmol, confId=cid)

    return rdkitmol, conformers


def plot_shapes(ratios, filename):
    """
    Plot molecule shapes.

    """

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(
        [ratios[i][1] for i in ratios],
        [ratios[i][0] for i in ratios],
        c='gold',
        edgecolors='white',
        alpha=1.0,
        s=80
    )

    ax.plot([0, 0.5, 1, 0], [1, 0.5, 1, 1], c='k', lw=2)
    ax.text(0.75, 1.03, 'sphere', fontsize=20)
    ax.text(0.4, 0.45, 'oblate', fontsize=20)
    ax.text(-0.05, 1.03, 'prolate', fontsize=20)

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('$I_1$ / $I_3$', fontsize=16)
    ax.set_ylabel('$I_2$ / $I_3$', fontsize=16)
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(0.4, 1.1)

    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')


def plot_diameters(diameters, filename):
    """
    Plot molecule diameters.

    """

    d1 = [diameters[i][0] for i in diameters]
    d2 = [diameters[i][1] for i in diameters]
    d3 = [diameters[i][2] for i in diameters]
    print(d1, d2, d3)

    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.4
    xlim = (0, round(max(d3)+1, 0))
    X_bins = np.arange(xlim[0], xlim[1], width)
    hist, bin_edges = np.histogram(a=d1, bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color='gold',
        edgecolor='k'
    )
    hist, bin_edges = np.histogram(a=d2, bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color='palegreen',
        edgecolor='k'
    )
    hist, bin_edges = np.histogram(a=d3, bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color='coral',
        edgecolor='k'
    )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'molecular diameter [$\mathrm{\AA}]', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')
