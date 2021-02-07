import pytest
import numpy as np
import molellipsize as mes


def test_get_inertial_ratios(molecule):
    mes_mol = mes.Molecule(molecule.molecule, molecule.conformers)
    conf_ratios = mes_mol.get_inertial_ratios()
    print(conf_ratios)
    test_ratio_1, test_ratio_2 = conf_ratios[0]
    assert np.isclose(
        test_ratio_1, molecule.inertial_ratios[0], rtol=0, atol=1E-1
    )
    assert np.isclose(
        test_ratio_2, molecule.inertial_ratios[1], rtol=0, atol=1E-1
    )


def test_get_ellipsoids(molecule):
    mes_mol = mes.Molecule(molecule.molecule, molecule.conformers)
    conf_ellipsoids = mes_mol.get_ellipsoids(1.0, 4.0, 0.5)

    test_center = conf_ellipsoids[0][0]
    for known, test_cent_dim in zip(
        molecule.ellipse_center,
        test_center
    ):
        assert np.isclose(known, test_cent_dim, rtol=0, atol=1E-1)

    test_diameters = conf_ellipsoids[0][1]
    print(test_diameters)
    for known, test_diam in zip(
        molecule.ellipse_diameters,
        test_diameters
    ):
        assert np.isclose(known, test_diam, rtol=0, atol=1E-5)
