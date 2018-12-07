# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import glob
import itertools
from os.path import join, splitext
import pytest
from pytest import approx
import numpy as np
import biotite
import biotite.structure.io.gro as gro
import biotite.structure.io.pdb as pdb
from biotite.structure.atoms import AtomArray
from .util import data_dir


@pytest.mark.parametrize(
    "path, single_model",
    itertools.product(
        glob.glob(join(data_dir, "*.gro")),
        [False, True]
    )
)
def test_array_conversion(path, single_model):
    model = 1 if single_model else None
    gro_file = gro.GROFile()
    gro_file.read(path)
    array1 = gro_file.get_structure(model=model)
    gro_file = gro.GROFile()
    gro_file.set_structure(array1)
    array2 = gro_file.get_structure(model=model)
    assert array1 == array2


@pytest.mark.parametrize("path", glob.glob(join(data_dir, "*.gro")))
def test_pdb_consistency(path):
    pdb_path = splitext(path)[0] + ".pdb"
    pdb_file = pdb.PDBFile()
    pdb_file.read(pdb_path)
    a1 = pdb_file.get_structure(model=1)
    gro_file = gro.GROFile()
    gro_file.read(path)
    a2 = gro_file.get_structure(model=1)

    assert a1.array_length() == a2.array_length()

    for category in ["res_id", "res_name", "atom_name"]:
        assert a1.get_annotation(category).tolist() == \
               a2.get_annotation(category).tolist()

    # Mind rounding errors when converting pdb to gro (A -> nm)
    assert a1.coord.flatten().tolist() \
        == approx(a2.coord.flatten().tolist(), abs=1e-2)


@pytest.mark.parametrize(
    "path, single_model",
    itertools.product(
        glob.glob(join(data_dir, "*.pdb")),
        [False, True]
    )
)
def test_pdb_to_gro(path, single_model):
    # Converting stacks between formats should not change data
    model = 1 if single_model else None
    
    # Read in data
    pdb_file = pdb.PDBFile()
    pdb_file.read(path)
    a1 = pdb_file.get_structure(model=model)

    # Save stack as gro
    tmp_file_name = biotite.temp_file("gro")
    gro_file = gro.GROFile()
    gro_file.set_structure(a1)
    gro_file.write(tmp_file_name)

    # Reload stack from gro
    gro_file = gro.GROFile()
    gro_file.read(tmp_file_name)
    a2 = gro_file.get_structure(model=model)

    assert a1.array_length() == a2.array_length()

    for category in ["res_id", "res_name", "atom_name"]:
        assert a1.get_annotation(category).tolist() == \
               a2.get_annotation(category).tolist()

    # Mind rounding errors when converting pdb to gro (A -> nm)
    assert a1.coord.flatten().tolist() \
        == approx(a2.coord.flatten().tolist(), abs=1e-2)


@pytest.mark.parametrize(
    "path, single_model",
    itertools.product(
        glob.glob(join(data_dir, "*.gro")),
        [False, True]
    )
)
def test_box_shape(path, single_model):
    model = 1 if single_model else None
    gro_file = gro.GROFile()
    gro_file.read(path)
    a = gro_file.get_structure(model=model)

    if isinstance(a, AtomArray):
        expected_box_dim = (3, 3)
    else:
        expected_box_dim = (len(a), 3, 3)
    assert expected_box_dim == a.box.shape


def test_box_parsing():
    path = join(data_dir, "1l2y.gro")
    gro_file = gro.GROFile()
    gro_file.read(path)
    a = gro_file.get_structure()
    expected_box = np.array([[
        [2.00000, 1.88562, 1.63299],
        [0.00000, 0.00000, 0.66667],
        [0.00000, -0.66667, 0.94281]
    ]])
    box = a.box
    assert np.array_equal(box, expected_box)



