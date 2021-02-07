"""
This module defines general-purpose objects, functions and classes.

"""

import numpy as np
from rdkit.Chem import AllChem as Chem
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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
    Plot molecule intertial ratios.

    Saves figure at `filename`.

    """

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(
        [ratios[i][0] for i in ratios],
        [ratios[i][1] for i in ratios],
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

    Saves figure at `filename`.

    """

    d1 = [diameters[i][0] for i in diameters]
    d2 = [diameters[i][1] for i in diameters]
    d3 = [diameters[i][2] for i in diameters]

    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.2
    xlim = (0, round(max(d3)+1, 0))
    X_bins = np.arange(xlim[0], xlim[1], width)
    hist, bin_edges = np.histogram(a=d1, bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color='none',
        edgecolor='gold',
        linewidth=3,
        label='d1',
    )
    hist, bin_edges = np.histogram(a=d2, bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color='none',
        edgecolor='palegreen',
        linewidth=3,
        label='d2',
    )
    hist, bin_edges = np.histogram(a=d3, bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color='none',
        edgecolor='coral',
        linewidth=3,
        label='d3',
    )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'molecular diameter [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(0, None)
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')


def plot_ellipsoid(
    center,
    diameter,
    rotation,
    plotAxes=False,
    cageColor='k',
    cageAlpha=0.2,
):
    """
    Plot an ellipsoid.

    Returns
    -------
    fig : class:`matplotlib.figure`
        Figure object with ellipsoid.

    ax : class:`matplotlib.axis`
        Axis object with ellipsoid.

    """

    radii = [i/2 for i in diameter]

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
        axes = np.array([
            [radii[0], 0.0, 0.0],
            [0.0, radii[1], 0.0],
            [0.0, 0.0, radii[2]]
        ])
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

    return fig, ax
