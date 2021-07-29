mol-ellipsize
=============

:author: Andrew Tarzia

Molecular size calculation based on ellipsoid fitting over N conformers.

This work was developed for the work in my repository `enzyme-screen` for screening enzymatic reactions for MOF@Enzyme systems.
Note that the implementation in this repository mirrors that in `enzyme-screen`, which also has the size calculation functionality, but is significantly more light-weight.

Please contact me with any questions (<andrew.tarzia@gmail.com>) or submit an issue!

.. image:: https://zenodo.org/badge/327267391.svg
   :target: https://zenodo.org/badge/latestdoi/327267391

Installation
------------

To get ``mol-ellipsize``, you can install it with pip::

    $ pip install mol-ellipsize

Make sure you also install rdkit, which is a dependency (version 2019.09.2.0 was used in `enzyme_screen`)::

    $ conda install -c rdkit rdkit=2020

Algorithm
---------

This code focusses on the calculation of the size of molecules within a conformer ensemble based on the fit of an ellipsoid around the molecules van der Waals cloud.
Any conformer ensemble can be provided through the rdkit .Molecule and Conformer classes.
However, helper functions are provided for generating ensembles using rdkit's ETKDG algorithm.

The ellipsoid fitting algorithm was modified from From: https://github.com/minillinim/ellipsoid.
The code is based on work by Nima Moshtagh <http://www.mathworks.com/matlabcentral/fileexchange/9542> and also by looking at <http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html>
It uses the Khachiyan algorithm to find the minimum volume ellipsoid that encompasses all points given to the function.

Examples
--------

The base example in ``examples/base_example.py`` shows the usage of this code to calculate the molecular size of molecules from SMILES strings.

A minimum example for calculating the size of 10 conformers of caffeine:

.. code-block:: python

    from rdkit.Chem import AllChem as Chem
    import molellipsize as mes


    rdkitmol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
    rdkitmol = Chem.AddHs(rdkitmol)
    Chem.SanitizeMol(rdkitmol)
    rdkitmol, conformers = mes.ETKDG_UFF_conformers(
        rdkitmol=rdkitmol,
        num_conformers=10,
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

*Adding coordinates and conformers without using RDKit*

The example in ``examples/arbitrary_coordinates.py`` shows the
fitting of an ellipsoid to arbitrary points with and without the
definition of a .Molecule.

Contributors and Acknowledgements
---------------------------------

I developed this code as a PhD student in the research groups of David Huang (<https://huanggroup.org/>) and Christian Doonan (<http://www.sumbydoonangroup.com/>) at the University of Adelaide.

License
-------

This project is licensed under the GPLv3 license.
