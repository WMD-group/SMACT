"""Test structure prediction module."""

from __future__ import annotations

import itertools
import json
import logging
import os
import pickle
from contextlib import contextmanager
from importlib.util import find_spec
from operator import itemgetter
from random import sample
from typing import TYPE_CHECKING, Any
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import requests
from pandas.testing import assert_frame_equal, assert_series_equal
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.structure_prediction.substitution_probability import (
    SubstitutionProbability,
)
from pymatgen.core import SETTINGS
from pymatgen.core import Lattice as PmgLattice
from pymatgen.core import Structure as PmgStructure

import smact.structure_prediction.structure as sp_struct
from smact import Species
from smact.structure_prediction import logger as sp_logger
from smact.structure_prediction.database import StructureDB, parse_mprest
from smact.structure_prediction.mutation import CationMutator
from smact.structure_prediction.prediction import StructurePredictor
from smact.structure_prediction.structure import SmactStructure

if TYPE_CHECKING:
    from collections.abc import Generator

MP_URL = "https://api.materialsproject.org"
MP_API_AVAILABLE = bool(find_spec("mp_api"))

try:
    skip_mprester_tests = requests.get(MP_URL, timeout=60).status_code != 200

except (ModuleNotFoundError, ImportError, requests.exceptions.RequestException):
    # Skip all MPRester tests if some downstream problem on the website, mp-api or whatever.
    skip_mprester_tests = True

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_STRUCT = os.path.join(files_dir, "test_struct")
TEST_POSCAR = os.path.join(files_dir, "test_poscar.txt")
TEST_PY_STRUCT = os.path.join(files_dir, "pymatgen_structure.json")
TEST_LAMBDA_JSON = os.path.join(files_dir, "test_lambda_tab.json")
TEST_LAMBDA_CSV = os.path.join(files_dir, "test_lambda_tab.csv")
TEST_PREDICTOR_DB = os.path.join(files_dir, "test_predictor.db")
TEST_PREDICTOR_TABLE = "TEST"

TEST_SPECIES: dict[str, list[tuple[str, int, int] | tuple[Species, int]]] = {
    "CaTiO3": [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)],
    "NaCl": [("Na", 1, 1), ("Cl", -1, 1)],
    "Fe": [("Fe", 0, 1)],
}


def generate_test_structure(comp: str) -> bool:
    """Generate a pickled test structure for comparison."""
    poscar_file = os.path.join(files_dir, f"{comp}.txt")
    s = SmactStructure.from_file(poscar_file)

    with open(TEST_STRUCT, "wb") as f:
        pickle.dump(s, f)

    with open(TEST_POSCAR, "w") as f:
        f.write(s.as_poscar())

    return True


@contextmanager
def ignore_warnings(logger: logging.Logger) -> Generator[int, None, None]:
    """Ignore logging warnings."""
    log_lvl_buff = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield log_lvl_buff
    logger.setLevel(log_lvl_buff)


def _assert_struct_almost_equal(s1: SmactStructure, s2: SmactStructure, places: int = 7):
    """
    Assert that two SmactStructures are almost equal.

    Almost equality dependent on how many decimal places the site coordinates
    are equal to.
    """
    must_equal = ["species", "lattice_param"]
    for cond in must_equal:
        assert getattr(s1, cond) == getattr(s2, cond)

    assert np.array_equal(s1.lattice_mat, s2.lattice_mat)
    assert list(s1.sites.keys()) == list(s2.sites.keys())

    for si, sj in zip(s1.sites.values(), s2.sites.values(), strict=True):
        for ci, cj in zip(si, sj, strict=True):
            for c1, c2 in zip(ci, cj, strict=True):
                assert c1 == pytest.approx(c2, abs=10**-places)


def _gen_empty_structure(species):
    """Generate an empty set of arguments for `SmactStructure` testing."""
    lattice_mat = np.array([[0] * 3] * 3)

    if isinstance(species[0][0], str):
        species_strs = [
            "{ele}{charge}{sign}".format(
                ele=spec[0],
                charge=abs(spec[1]),
                sign="+" if spec[1] >= 0 else "-",
            )
            for spec in species
        ]
    else:
        species_strs = [
            "{ele}{charge}{sign}".format(
                ele=spec[0].symbol,
                charge=abs(spec[0].oxidation),
                sign="+" if spec[0].oxidation >= 0 else "-",
            )
            for spec in species
        ]

    sites = {spec: [[]] for spec in species_strs}
    return species, lattice_mat, sites


# ---------------------------------------------------------------------------
# StructureTest
# ---------------------------------------------------------------------------


def test_as_poscar():
    for comp in TEST_SPECIES:
        comp_file = os.path.join(files_dir, f"{comp}.txt")
        with open(comp_file) as f:
            struct = SmactStructure.from_file(comp_file)
            assert struct.as_poscar() == f.read()


def test_from_py_struct():
    with open(TEST_PY_STRUCT) as f:
        d = json.load(f)
        py_structure = PmgStructure.from_dict(d)

    with ignore_warnings(sp_logger):
        s1 = SmactStructure.from_py_struct(py_structure)

    s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))

    _assert_struct_almost_equal(s1, s2)


def test_from_py_struct_icsd():
    with open(TEST_PY_STRUCT) as f:
        d = json.load(f)
        py_structure = PmgStructure.from_dict(d)

    with ignore_warnings(sp_logger):
        s1 = SmactStructure.from_py_struct(py_structure, determine_oxi="comp_ICSD")

    s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))

    _assert_struct_almost_equal(s1, s2)


def test_has_species():
    s1 = SmactStructure(*_gen_empty_structure([("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]))

    assert s1.has_species(("Ba", 2))
    assert not s1.has_species(("Ba", 3))
    assert not s1.has_species(("Ca", 2))


def test_smactStruc_comp_key():
    s1 = SmactStructure(*_gen_empty_structure([("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]))
    s2 = SmactStructure(*_gen_empty_structure([("Fe", 2, 1), ("Fe", 3, 2), ("O", -2, 4)]))

    Ba = Species("Ba", 2)
    ox = Species("O", -2)
    F = Species("F", -1)
    Fe2 = Species("Fe", 2)
    Fe3 = Species("Fe", 3)

    s3 = SmactStructure(*_gen_empty_structure([(Ba, 2), (ox, 1), (F, 2)]))
    s4 = SmactStructure(*_gen_empty_structure([(Fe2, 1), (Fe3, 2), (ox, 4)]))

    Ba_2OF_2 = "Ba_2_2+F_2_1-O_1_2-"
    Fe_3O_4 = "Fe_2_3+Fe_1_2+O_4_2-"
    assert s1.composition() == Ba_2OF_2
    assert s2.composition() == Fe_3O_4
    assert s3.composition() == Ba_2OF_2
    assert s4.composition() == Fe_3O_4


def test_smactStruc_from_file():
    with open(TEST_STRUCT, "rb") as f:
        s1 = pickle.load(f)

    s2 = SmactStructure.from_file(TEST_POSCAR)

    assert s1 == s2


def test_equality():
    struct_files = [os.path.join(files_dir, f"{x}.txt") for x in ["CaTiO3", "NaCl"]]
    CaTiO3 = SmactStructure.from_file(struct_files[0])
    NaCl = SmactStructure.from_file(struct_files[1])

    same_ref = CaTiO3
    assert CaTiO3 == same_ref
    assert CaTiO3 != "CaTiO3"
    assert CaTiO3 != NaCl


def test_ele_stoics():
    s1 = SmactStructure(*_gen_empty_structure([("Fe", 2, 1), ("Fe", 3, 2), ("O", -2, 4)]))
    s1_stoics = {"Fe": 3, "O": 4}
    s2 = SmactStructure(*_gen_empty_structure([("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]))
    s2_stoics = {"Ba": 2, "O": 1, "F": 2}

    for test, expected in [(s1, s1_stoics), (s2, s2_stoics)]:
        assert SmactStructure._get_ele_stoics(test.species) == expected


@pytest.mark.integration
@pytest.mark.skipif(
    (skip_mprester_tests or not MP_API_AVAILABLE or not (os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY"))),
    reason="Materials Project API not available or not configured.",
)
def test_from_mp():
    api_key = os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY")

    for comp, species in TEST_SPECIES.items():
        comp_file = os.path.join(files_dir, f"{comp}.txt")
        local_struct = SmactStructure.from_file(comp_file)
        mp_struct = SmactStructure.from_mp(species, api_key)
        assert local_struct == mp_struct


def test_hash():
    s1 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
    s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
    assert isinstance(hash(s1), int)
    assert hash(s1) == hash(s2)
    # Different structures must produce different hashes
    s3 = SmactStructure.from_file(os.path.join(files_dir, "NaCl.txt"))
    assert hash(s1) != hash(s3)


def test_sanitise_species_error_branches():
    bad_not_a_list: Any = "not_a_list"
    with pytest.raises(TypeError):  # line 199: not a list
        SmactStructure._sanitise_species(bad_not_a_list)
    with pytest.raises(ValueError):  # line 201: empty list
        SmactStructure._sanitise_species([])
    bad_non_tuples: Any = ["not_a_tuple"]
    with pytest.raises(TypeError):  # line 203: list of non-tuples
        SmactStructure._sanitise_species(bad_non_tuples)
    bad_wrong_length: Any = [("a", "b", "c", "d")]
    with pytest.raises(ValueError):  # line 211: tuple has wrong length
        SmactStructure._sanitise_species(bad_wrong_length)
    bad_first_element: Any = [(42, 2, 1)]
    with pytest.raises(TypeError):  # line 223: first element is not str or Species
        SmactStructure._sanitise_species(bad_first_element)


def test_parse_py_sites_type_error():
    bad_input: Any = "not_a_structure"
    _parse_py_sites = "_SmactStructure__parse_py_sites"
    with pytest.raises(TypeError, match=r"pymatgen\.core\.Structure"):
        getattr(SmactStructure, _parse_py_sites)(bad_input)


def test_parse_py_sites_unit_charge():
    # NaCl decorated by BV gives Na+ and Cl- — unit charges trigger line 262
    nacl = PmgStructure.from_spacegroup(
        "Fm-3m",
        PmgLattice.cubic(5.6),
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    with ignore_warnings(sp_logger):
        s = SmactStructure.from_py_struct(nacl, determine_oxi="BV")
    assert isinstance(s, SmactStructure)
    # Confirm Na+1 and Cl-1 are present
    species_list = [sp[0] for sp in s.species]
    assert "Na" in species_list
    assert "Cl" in species_list


def test_parse_py_sites_zero_charge():
    # A plain (undecorated) Structure used with determine_oxi='predecorated'
    # yields species_string = 'Fe' (no digits) → else branch
    fe_struct = PmgStructure(PmgLattice.cubic(2.87), ["Fe"], [[0, 0, 0]])
    with ignore_warnings(sp_logger):
        s = SmactStructure.from_py_struct(fe_struct, determine_oxi="predecorated")
    assert isinstance(s, SmactStructure)
    assert s.species[0][1] == 0  # charge = 0


def test_from_py_struct_type_error():
    bad_input: Any = "not_a_structure"
    with pytest.raises(TypeError, match=r"pymatgen\.core\.Structure"):
        SmactStructure.from_py_struct(bad_input)


def test_from_py_struct_both_success():
    with open(TEST_PY_STRUCT) as f:
        py_structure = PmgStructure.from_dict(json.load(f))
    with ignore_warnings(sp_logger):
        s = SmactStructure.from_py_struct(py_structure, determine_oxi="both")
    assert isinstance(s, SmactStructure)


def test_from_py_struct_both_bv_fallback():
    with open(TEST_PY_STRUCT) as f:
        py_structure = PmgStructure.from_dict(json.load(f))

    with (
        patch.object(BVAnalyzer, "get_oxi_state_decorated_structure", side_effect=ValueError("BV failed")),
        ignore_warnings(sp_logger),
    ):
        s = SmactStructure.from_py_struct(py_structure, determine_oxi="both")
    assert isinstance(s, SmactStructure)


def test_from_py_struct_predecorated():
    with open(TEST_PY_STRUCT) as f:
        py_structure = PmgStructure.from_dict(json.load(f))
    decorated = BVAnalyzer().get_oxi_state_decorated_structure(py_structure)
    with ignore_warnings(sp_logger):
        s = SmactStructure.from_py_struct(decorated, determine_oxi="predecorated")
    assert isinstance(s, SmactStructure)


def test_from_py_struct_invalid_determine_oxi():
    with open(TEST_PY_STRUCT) as f:
        py_structure = PmgStructure.from_dict(json.load(f))
    with pytest.raises(ValueError, match="determine_oxi"):
        SmactStructure.from_py_struct(py_structure, determine_oxi="unknown")


def test_from_mp_no_api_key():
    env_keys = ["MP_API_KEY", "PMG_MAPI_KEY"]
    saved = {k: os.environ.pop(k, None) for k in env_keys}
    try:
        # Replace SETTINGS with an empty dict so .get("PMG_MAPI_KEY") returns None
        with (
            patch.object(sp_struct, "SETTINGS", {}),
            pytest.raises(ValueError, match="No Materials Project API key"),
        ):
            SmactStructure.from_mp([("Na", 1, 1), ("Cl", -1, 1)])
    finally:
        for k, v in saved.items():
            if v is not None:
                os.environ[k] = v


def test_from_mp_mocked_new_api():
    with open(TEST_PY_STRUCT) as f:
        real_struct = PmgStructure.from_dict(json.load(f))

    mock_cm = MagicMock()
    mock_cm.materials.summary.search.return_value = [real_struct]
    mock_rester = MagicMock()
    mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
    mock_rester.return_value.__exit__ = MagicMock(return_value=False)

    species: list[tuple[str, int, int] | tuple[Species, int]] = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
    with (
        patch.object(sp_struct, "MPResterNew", new=mock_rester),
        patch.object(sp_struct, "HAS_MP_API", new=True),
        ignore_warnings(sp_logger),
    ):
        s = SmactStructure.from_mp(species, api_key="x" * 32)
    assert isinstance(s, SmactStructure)


def test_from_mp_mocked_legacy_api():
    with open(TEST_PY_STRUCT) as f:
        real_struct = PmgStructure.from_dict(json.load(f))

    mock_cm = MagicMock()
    mock_cm.query.return_value = [{"structure": real_struct}]
    mock_rester = MagicMock()
    mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
    mock_rester.return_value.__exit__ = MagicMock(return_value=False)

    species: list[tuple[str, int, int] | tuple[Species, int]] = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
    with (
        patch.object(sp_struct, "MPRester", new=mock_rester),
        patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=True),
        ignore_warnings(sp_logger),
    ):
        s = SmactStructure.from_mp(species, api_key="short_key")
    assert isinstance(s, SmactStructure)


def test_from_mp_mocked_legacy_32char():
    with open(TEST_PY_STRUCT) as f:
        real_struct = PmgStructure.from_dict(json.load(f))

    mock_cm = MagicMock()
    mock_cm.get_structures.return_value = [real_struct]
    mock_rester = MagicMock()
    mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
    mock_rester.return_value.__exit__ = MagicMock(return_value=False)

    species: list[tuple[str, int, int] | tuple[Species, int]] = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
    with (
        patch.object(sp_struct, "MPRester", new=mock_rester),
        patch.object(sp_struct, "HAS_MP_API", new=False),
        patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=True),
        ignore_warnings(sp_logger),
    ):
        s = SmactStructure.from_mp(species, api_key="x" * 32)
    assert isinstance(s, SmactStructure)


def test_from_mp_determine_oxi_branches():
    with open(TEST_PY_STRUCT) as f:
        real_struct = PmgStructure.from_dict(json.load(f))

    def make_mock_rester(return_val):
        mock_cm = MagicMock()
        mock_cm.materials.summary.search.return_value = [return_val]
        mock_rester = MagicMock()
        mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
        mock_rester.return_value.__exit__ = MagicMock(return_value=False)
        return mock_rester

    species: list[tuple[str, int, int] | tuple[Species, int]] = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
    api_key_32 = "x" * 32

    # comp_ICSD branch (lines 420-424)
    with (
        patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
        patch.object(sp_struct, "HAS_MP_API", new=True),
        ignore_warnings(sp_logger),
    ):
        s = SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="comp_ICSD")
        assert isinstance(s, SmactStructure)

    # both branch — BV succeeds (lines 426-430)
    with (
        patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
        patch.object(sp_struct, "HAS_MP_API", new=True),
        ignore_warnings(sp_logger),
    ):
        s = SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="both")
        assert isinstance(s, SmactStructure)

    # both branch — BV fails, falls back to comp_ICSD (lines 431-435)
    with (
        patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
        patch.object(sp_struct, "HAS_MP_API", new=True),
        patch.object(BVAnalyzer, "get_oxi_state_decorated_structure", side_effect=ValueError("BV failed")),
        ignore_warnings(sp_logger),
    ):
        s = SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="both")
        assert isinstance(s, SmactStructure)

    # invalid determine_oxi → ValueError (lines 436-439)
    with (
        patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
        patch.object(sp_struct, "HAS_MP_API", new=True),
        pytest.raises(ValueError, match="determine_oxi"),
    ):
        SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="invalid_oxi")


def test_from_mp_empty_result_raises():
    mock_cm = MagicMock()
    mock_cm.materials.summary.search.return_value = []
    mock_rester = MagicMock()
    mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
    mock_rester.return_value.__exit__ = MagicMock(return_value=False)

    with (
        patch.object(sp_struct, "MPResterNew", new=mock_rester),
        patch.object(sp_struct, "HAS_MP_API", new=True),
        pytest.raises(ValueError, match="Could not find composition"),
    ):
        SmactStructure.from_mp([("Na", 1, 1), ("Cl", -1, 1)], api_key="x" * 32)


def test_from_mp_no_legacy_mprester_short_key():
    species: list[tuple[str, int, int] | tuple[Species, int]] = [("Na", 1, 1), ("Cl", -1, 1)]
    with (
        patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=False),
        pytest.raises(ImportError, match="pymatgen legacy MPRester is not available"),
    ):
        SmactStructure.from_mp(species, api_key="short_key")


def test_from_mp_no_legacy_mprester_no_mp_api():
    species: list[tuple[str, int, int] | tuple[Species, int]] = [("Na", 1, 1), ("Cl", -1, 1)]
    api_key_32 = "x" * 32
    with (
        patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=False),
        patch.object(sp_struct, "HAS_MP_API", new=False),
        pytest.raises(ImportError, match="Neither mp-api nor pymatgen legacy MPRester"),
    ):
        SmactStructure.from_mp(species, api_key=api_key_32)


def test_as_py_struct():
    s1 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
    s1_pym = s1.as_py_struct()
    with open(TEST_PY_STRUCT) as f:
        d = json.load(f)
    py_structure = PmgStructure.from_dict(d)

    assert s1_pym == py_structure


def test_reduced_formula():
    s1 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
    s2 = SmactStructure.from_file(os.path.join(files_dir, "NaCl.txt"))

    assert s1.reduced_formula() == "CaTiO3"
    assert s2.reduced_formula() == "NaCl"


# ---------------------------------------------------------------------------
# StructureDBTest
# ---------------------------------------------------------------------------


def test_db_interface(tmp_path):
    test_db = str(tmp_path / "test_db.tmp")
    test_table = "Structures"
    test_mp_table = "Structures1"

    db = StructureDB(test_db)

    db.add_table(test_table)

    struct_file = os.path.join(files_dir, "CaTiO3.txt")
    struct = SmactStructure.from_file(struct_file)

    db.add_struct(struct, test_table)

    struct_list = db.get_structs(struct.composition(), test_table)
    assert len(struct_list) == 1
    assert struct_list[0] == struct

    struct_files = [os.path.join(files_dir, f"{x}.txt") for x in ["NaCl", "Fe"]]
    structs = [SmactStructure.from_file(fname) for fname in struct_files]

    db.add_structs(structs, test_table)

    test_with_species_args = [
        [("Na", 1)],
        [("Cl", -1)],
        [("Na", 1), ("Cl", -1)],
        [("Cl", -1), ("Na", 1)],
        [("Cl", -1)],
        [("Na", 1), ("Cl", 1)],
        [("O", -2)],
        [("Ca", 2), ("Ti", 4), ("O", -2)],
    ]

    test_with_species_exp = [
        [structs[0]],
        [structs[0]],
        [structs[0]],
        [structs[0]],
        [structs[0]],
        [],
        [struct],
        [struct],
    ]

    for spec, expected in zip(test_with_species_args, test_with_species_exp, strict=False):
        assert db.get_with_species(spec, test_table) == expected

    mp_strucs = [PmgStructure.from_file(os.path.join(files_dir, f)) for f in ["CaTiO3.json", "NaCl.json", "Fe3O4.json"]]
    mp_data = [
        {"material_id": mpid, "structure": s}
        for mpid, s in zip(["mp-4019", "mp-22862", "mp-19306"], mp_strucs, strict=False)
    ]

    added: int = db.add_mp_icsd(test_mp_table, mp_data)
    assert added == 3


def test_db_rollback_on_exception(tmp_path):
    def _raise_intentional_error():
        msg = "intentional error for rollback"
        raise RuntimeError(msg)

    db_path = str(tmp_path / "rollback.db")

    db = StructureDB(db_path)
    db.add_table("rollback_test")

    try:
        with db as c:
            c.execute("INSERT INTO rollback_test VALUES (?, ?)", ("comp", "poscar"))
            _raise_intentional_error()
    except RuntimeError:
        pass

    # Rollback means the INSERT should not have been committed
    with db as c:
        c.execute("SELECT * FROM rollback_test")
        rows = c.fetchall()
    assert len(rows) == 0


def test_add_structs_with_none(tmp_path):
    db_path = str(tmp_path / "test_none.db")
    db = StructureDB(db_path)
    db.add_table("test_none")
    struct = SmactStructure.from_file(os.path.join(files_dir, "NaCl.txt"))
    added = db.add_structs([None, struct, None], "test_none")
    assert added == 1


def test_parse_mprest_exception_path():
    # Passing a non-Structure as 'structure' causes from_py_struct → TypeError → caught
    parse_mprest({"structure": "not_a_structure", "material_id": "mp-test"})


# ---------------------------------------------------------------------------
# CationMutatorTest
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cation_mutator_data():
    test_struct = SmactStructure.from_file(TEST_POSCAR)

    test_mutator = CationMutator.from_json(lambda_json=TEST_LAMBDA_JSON)
    test_pymatgen_mutator = CationMutator.from_json(lambda_json=None, alpha=lambda _x, _y: -5)

    # 5 random test species -> 5! test pairs
    test_species = sample(list(test_pymatgen_mutator.specs), 5)
    test_pairs = list(itertools.combinations_with_replacement(test_species, 2))

    pymatgen_sp = SubstitutionProbability(lambda_table=None, alpha=-5)

    return {
        "test_struct": test_struct,
        "test_mutator": test_mutator,
        "test_pymatgen_mutator": test_pymatgen_mutator,
        "test_species": test_species,
        "test_pairs": test_pairs,
        "pymatgen_sp": pymatgen_sp,
    }


def test_lambda_tab_pop(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    lambda_dat = [
        [-5.0, 0.5, -5.0],
        [0.5, -5.0, 0.3],
        [-5.0, 0.3, -5.0],
    ]
    labels = ["A", "B", "C"]

    exp_lambda = pd.DataFrame(lambda_dat, index=pd.Index(labels), columns=pd.Index(labels))

    assert_frame_equal(
        test_mutator.lambda_tab,
        exp_lambda,
        check_names=False,
    )


def test_partition_func_Z(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    # 2e^0.5 + 2e^0.3 + 5e^{-5} \approx 6.0308499
    assert pytest.approx(6.0308499) == test_mutator.Z


def test_pymatgen_lambda_import(cation_mutator_data):
    test_pymatgen_mutator = cation_mutator_data["test_pymatgen_mutator"]
    assert isinstance(test_pymatgen_mutator.lambda_tab, pd.DataFrame)


def test_lambda_interface(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    test_cases = [itertools.permutations(x) for x in [("A", "B"), ("A", "C"), ("B", "C")]]

    expected = [0.5, -5.0, 0.3]

    for test_case, expectation in zip(test_cases, expected, strict=False):
        for spec_comb in test_case:
            s1, s2 = spec_comb
            assert test_mutator.get_lambda(s1, s2) == expectation


def test_ion_mutation(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    ca_file = os.path.join(files_dir, "CaTiO3.txt")
    ba_file = os.path.join(files_dir, "BaTiO3.txt")

    CaTiO3 = SmactStructure.from_file(ca_file)
    BaTiO3 = SmactStructure.from_file(ba_file)

    mutation = test_mutator._mutate_structure(CaTiO3, "Ca2+", "Ba2+")
    assert mutation == BaTiO3

    na_file = os.path.join(files_dir, "NaCl.txt")
    NaCl = SmactStructure.from_file(na_file)

    with pytest.raises(ValueError):
        test_mutator._mutate_structure(NaCl, "Na1+", "Na2+")


def test_sub_prob(cation_mutator_data):
    test_pymatgen_mutator = cation_mutator_data["test_pymatgen_mutator"]
    test_pairs = cation_mutator_data["test_pairs"]
    pymatgen_sp = cation_mutator_data["pymatgen_sp"]
    for s1, s2 in test_pairs:
        assert pymatgen_sp.prob(s1, s2) == pytest.approx(test_pymatgen_mutator.sub_prob(s1, s2))


def test_cond_sub_probs(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    for s1 in ["A", "B", "C"]:
        cond_sub_probs_test = test_mutator.cond_sub_probs(s1)

        vals = [(s1, s2, test_mutator.cond_sub_prob(s1, s2)) for s2 in ["A", "B", "C"]]

        test_df = pd.DataFrame(vals)
        test_df: pd.DataFrame = test_df.pivot_table(index=0, columns=1, values=2)

        # Slice to convert to series
        assert_series_equal(cond_sub_probs_test, test_df.iloc[0])


def test_cond_sub_prob(cation_mutator_data):
    test_pymatgen_mutator = cation_mutator_data["test_pymatgen_mutator"]
    test_pairs = cation_mutator_data["test_pairs"]
    pymatgen_sp = cation_mutator_data["pymatgen_sp"]
    for s1, s2 in test_pairs:
        assert pymatgen_sp.cond_prob(s1, s2) == pytest.approx(test_pymatgen_mutator.cond_sub_prob(s1, s2))


def test_pair_corr(cation_mutator_data):
    test_pymatgen_mutator = cation_mutator_data["test_pymatgen_mutator"]
    test_pairs = cation_mutator_data["test_pairs"]
    pymatgen_sp = cation_mutator_data["pymatgen_sp"]
    for s1, s2 in test_pairs:
        assert pymatgen_sp.pair_corr(s1, s2) == pytest.approx(test_pymatgen_mutator.pair_corr(s1, s2))


def test_from_df(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    lambda_df = pd.read_csv(TEST_LAMBDA_CSV, index_col=0)

    csv_test = CationMutator(lambda_df=lambda_df)

    assert_frame_equal(
        csv_test.lambda_tab,
        test_mutator.lambda_tab,
        check_names=False,
    )


def test_complete_cond_probs(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    pairs = itertools.product(["A", "B", "C"], repeat=2)

    vals = [
        (
            s1,
            s2,
            test_mutator.cond_sub_prob(s1, s2),
        )
        for s1, s2 in pairs
    ]

    cond_probs = pd.DataFrame(vals)
    cond_probs = cond_probs.pivot_table(index=0, columns=1, values=2)

    assert_frame_equal(test_mutator.complete_cond_probs(), cond_probs)


def test_complete_sub_probs(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    pairs = itertools.product(["A", "B", "C"], repeat=2)

    vals = [
        (
            s1,
            s2,
            test_mutator.sub_prob(s1, s2),
        )
        for s1, s2 in pairs
    ]

    sub_probs = pd.DataFrame(vals)
    sub_probs = sub_probs.pivot_table(index=0, columns=1, values=2)

    assert_frame_equal(test_mutator.complete_sub_probs(), sub_probs)


def test_complete_pair_corrs(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    pairs = itertools.product(["A", "B", "C"], repeat=2)

    vals = [
        (
            s1,
            s2,
            test_mutator.pair_corr(s1, s2),
        )
        for s1, s2 in pairs
    ]

    pair_corrs = pd.DataFrame(vals)
    pair_corrs = pair_corrs.pivot_table(index=0, columns=1, values=2)

    assert_frame_equal(test_mutator.complete_pair_corrs(), pair_corrs)


def test_same_spec_probs(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    probs = test_mutator.same_spec_probs()
    assert isinstance(probs, pd.Series)
    assert len(probs) == 3
    # All values should be positive
    assert (probs > 0).all()


def test_same_spec_cond_probs(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    probs = test_mutator.same_spec_cond_probs()
    # Conditional probabilities should be between 0 and 1
    assert (probs >= 0).all()
    assert (probs <= 1).all()


def test_get_lambda_unknown_species(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    lam = test_mutator.get_lambda("X", "Y")
    assert lam == -5.0


def test_get_lambdas_unknown_species_raises(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    with pytest.raises(ValueError, match="not in lambda table"):
        test_mutator.get_lambdas("UNKNOWN_SPEC")


def test_populate_lambda_alpha_none_raises():
    lambda_df = pd.DataFrame(
        [[1.0, np.nan], [np.nan, 1.0]],
        index=pd.Index(["A", "B"]),
        columns=pd.Index(["A", "B"]),
    )
    with pytest.raises(ValueError, match="alpha function must not be None"):
        CationMutator(lambda_df, alpha=None)


def test_populate_lambda_asymmetric_mirror():
    # Build a DataFrame where "A" appears only as a row index and "B"
    # only as a column, so (A, B) exists but (B, A) raises KeyError.
    # A third species "C" keeps the table from being trivially complete.
    lambda_df = pd.DataFrame(
        {"B": {"A": 2.0, "C": np.nan}, "C": {"A": np.nan, "C": -1.0}},
    )
    # (B, A) is KeyError, (A, B) = 2.0 → mirror_lambda(B, A) on line 147
    # (B, C) is KeyError, (C, B) is also KeyError → add_alpha on line 155
    # But we need the *outer* KeyError branch to mirror.
    # Construct so (s1, s2) = ("B", "A") raises KeyError → outer except,
    # then (s2, s1) = ("A", "B") = 2.0 → mirror_lambda (lines 152-153).
    cm = CationMutator(lambda_df, alpha=lambda s1, s2: -5.0)

    # After population, (B, A) should be mirrored from (A, B)
    assert cm.lambda_tab.loc["B", "A"] == 2.0
    assert cm.lambda_tab.loc["A", "B"] == 2.0


def test_nary_mutate_structure():
    ca_file = os.path.join(files_dir, "CaTiO3.txt")
    CaTiO3 = SmactStructure.from_file(ca_file)

    # Replace Ca2+ → Ba2+ and Ti4+ → Zr4+ (both same-charge substitutions)
    mutated = CationMutator._nary_mutate_structure(CaTiO3, ["Ca2+", "Ti4+"], ["Ba2+", "Zr4+"])
    assert isinstance(mutated, SmactStructure)
    spec_strs = mutated.get_spec_strs()
    assert "Ba2+" in spec_strs
    assert "Zr4+" in spec_strs
    assert "Ca2+" not in spec_strs
    assert "Ti4+" not in spec_strs


def test_nary_mutate_structure_not_neutral():
    ca_file = os.path.join(files_dir, "CaTiO3.txt")
    CaTiO3 = SmactStructure.from_file(ca_file)

    # Replace Ca2+ → Na1+ breaks charge neutrality
    with pytest.raises(ValueError, match="not charge neutral"):
        CationMutator._nary_mutate_structure(CaTiO3, ["Ca2+"], ["Na1+"])


def test_unary_substitute(cation_mutator_data):
    test_pymatgen_mutator = cation_mutator_data["test_pymatgen_mutator"]
    ca_file = os.path.join(files_dir, "CaTiO3.txt")
    CaTiO3 = SmactStructure.from_file(ca_file)
    # Use pymatgen mutator which has a full lambda table
    results = list(test_pymatgen_mutator.unary_substitute(CaTiO3, thresh=1e-3))
    # Should yield at least some substitutions
    assert len(results) > 0
    for struct, prob, old_spec, new_spec in results:
        assert isinstance(struct, SmactStructure)
        assert prob > 0
        assert old_spec != new_spec


def test_sub_probs(cation_mutator_data):
    test_mutator = cation_mutator_data["test_mutator"]
    probs = test_mutator.sub_probs("A")
    assert isinstance(probs, pd.Series)
    assert (probs >= 0).all()


# ---------------------------------------------------------------------------
# PredictorTest
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def predictor_data():
    db = StructureDB(TEST_PREDICTOR_DB)
    # NOTE: This may break if the pymatgen lambda table is updated
    cm = CationMutator.from_json()
    table = TEST_PREDICTOR_TABLE
    return {"db": db, "cm": cm, "table": table}


def test_prediction(predictor_data):
    db = predictor_data["db"]
    cm = predictor_data["cm"]
    table = predictor_data["table"]

    sp = StructurePredictor(cm, db, table)

    test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]

    # At time of creation, results with thresh=0 and include_same=True are as follows:
    #
    # Ca_1_2+O_4_2-Ti_2_4+: p=1.0, from Ca_1_2+O_4_2-Ti_2_4+
    # Ca_1_2+O_3_2-Ti_1_4+: p=0.0088406421433475, from Ca_1_2+O_3_2-Si_1_4+
    # Ca_1_2+O_3_2-Ti_1_4+: p=0.05502122697602804, from Ca_1_2+Mo_1_4+O_3_2-
    # Ca_1_2+O_3_2-Ti_1_4+: p=0.0445515892332118, from Ca_1_2+O_3_2-V_1_4+
    # Ca_2_2+O_4_2-Ti_1_4+: p=0.0274412135813993, from Fe_2_2+O_4_2-Ti_1_4+

    predictions = list(sp.predict_structs(test_specs, thresh=0.02, include_same=False))

    expected_comps = [
        "Ca_2_2+O_8_2-Sn_3_4+",
        "Ca_1_2+Mo_1_4+O_3_2-",
        "Ca_1_2+O_3_2-V_1_4+",
        "Ca_2_2+O_8_2-V_3_4+",
        "Fe_2_2+O_4_2-Ti_1_4+",
    ]

    expected_probs = [
        0.05604346851066936,
        0.05502122697602808,
        0.04455158923321182,
        0.04455158923321182,
        0.02744121358139931,
    ]

    # Sort from highest to lowest probability
    predictions.sort(key=itemgetter(1), reverse=True)

    parents = [x[2].composition() for x in predictions]

    assert parents == expected_comps

    probs = [x[1] for x in predictions]

    for p1, p2 in zip(probs, expected_probs, strict=False):
        assert p1 == pytest.approx(p2)


def test_prediction_include_same(predictor_data):
    db = predictor_data["db"]
    cm = predictor_data["cm"]
    table = predictor_data["table"]

    sp = StructurePredictor(cm, db, table)
    test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
    predictions = list(sp.predict_structs(test_specs, thresh=0.0, include_same=True))
    # The first results should be identical structures with probability 1.0
    same_structs = [(s, p, parent) for s, p, parent in predictions if p == 1.0]
    assert len(same_structs) > 0
    for s, _p, parent in same_structs:
        assert s == parent


def test_prediction_thresh_none(predictor_data):
    db = predictor_data["db"]
    cm = predictor_data["cm"]
    table = predictor_data["table"]

    sp = StructurePredictor(cm, db, table)
    test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
    predictions = list(sp.predict_structs(test_specs, thresh=None, include_same=False))
    # With no threshold, we should get at least as many as with a threshold
    predictions_thresh = list(sp.predict_structs(test_specs, thresh=0.02, include_same=False))
    assert len(predictions) >= len(predictions_thresh)


def test_nary_predict_n_ary_none(predictor_data):
    db = predictor_data["db"]
    cm = predictor_data["cm"]
    table = predictor_data["table"]

    sp = StructurePredictor(cm, db, table)
    test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
    predictions = list(sp.nary_predict_structs(test_specs, n_ary=None, include_same=True))
    # Only yields identical structures (n_ary=None → early return)
    for _s, p, _parent in predictions:
        assert p == 1.0


def test_nary_predict_equal_species_returns_early(predictor_data):
    db = predictor_data["db"]
    cm = predictor_data["cm"]
    table = predictor_data["table"]

    sp = StructurePredictor(cm, db, table)
    test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
    # n_ary=3 means we'd replace all 3 species, leaving 0 common → early return
    predictions = list(sp.nary_predict_structs(test_specs, n_ary=3, include_same=False))
    assert len(predictions) == 0


def test_nary_predict_include_same_false(predictor_data):
    db = predictor_data["db"]
    cm = predictor_data["cm"]
    table = predictor_data["table"]

    sp = StructurePredictor(cm, db, table)
    test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
    predictions = list(sp.nary_predict_structs(test_specs, n_ary=2, include_same=False))
    for _s, p, _parent in predictions:
        assert p != 1.0


# ---------------------------------------------------------------------------
# PredictorCoverageTest
# ---------------------------------------------------------------------------


def _make_structure(species_tuples):
    """Create a minimal SmactStructure for DB storage.

    species_tuples: list of (element, charge, stoichiometry).
    Each species gets stoichiometry-many sites at arbitrary coordinates.
    """
    lattice_mat = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=float)
    sites = {}
    for i, (ele, charge, stoic) in enumerate(species_tuples):
        sign = "+" if charge >= 0 else "-"
        spec_str = f"{ele}{abs(charge)}{sign}"
        sites[spec_str] = [[float(j), float(i), 0.0] for j in range(stoic)]
    return SmactStructure(species_tuples, lattice_mat, sites)


@pytest.fixture(scope="module")
def predictor_coverage_data(tmp_path_factory):
    table = "COV"
    tmp_path = tmp_path_factory.mktemp("predictor_coverage")
    db_path = str(tmp_path / "coverage.db")

    nacl = _make_structure([("Na", 1, 1), ("Cl", -1, 1)])
    ca_ti_o3 = _make_structure([("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)])
    # Parent with Fe4+ — NOT in the lambda table → triggers KeyError path
    ca_fe_o3 = _make_structure([("Ca", 2, 1), ("Fe", 4, 1), ("O", -2, 3)])
    # 4-species structure for n_ary=3 test
    ba_ca_ti_o4 = _make_structure([("Ba", 2, 1), ("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 4)])
    # Non-neutral structure: Na1+(x2) Cl1-(x1) → charge +1
    non_neutral = _make_structure([("Na", 1, 2), ("Cl", -1, 1)])

    db = StructureDB(db_path)
    db.add_table(table)
    for s in [nacl, ca_ti_o3, ca_fe_o3, ba_ca_ti_o4, non_neutral]:
        db.add_struct(s, table)

    # Minimal CationMutator — Fe4+ intentionally excluded
    specs = ["Na1+", "Cl1-", "Ca2+", "Ti4+", "O2-", "Ba2+", "Si4+", "K1+"]
    n = len(specs)
    data = np.full((n, n), -5.0)
    np.fill_diagonal(data, 0.0)
    lambda_df = pd.DataFrame(data, index=pd.Index(specs), columns=pd.Index(specs))
    cm = CationMutator(lambda_df)

    return {"db": db, "cm": cm, "table": table}


def test_predict_duplicate_species_continue(predictor_coverage_data):
    db = predictor_coverage_data["db"]
    cm = predictor_coverage_data["cm"]
    table = predictor_coverage_data["table"]

    sp = StructurePredictor(cm, db, table)
    species = [("Na", 1), ("Na", 1), ("Cl", -1)]
    predictions = list(sp.predict_structs(species, thresh=0.0, include_same=False))
    assert isinstance(predictions, list)


def test_predict_parent_multi_extra_species(predictor_coverage_data):
    cm = predictor_coverage_data["cm"]
    table = predictor_coverage_data["table"]

    species = [("Ca", 2), ("Ti", 4), ("O", -2)]
    # Parent with Ca2+ (matching) + Fe2+, Mn4+ (both differ from target)
    parent = _make_structure([("Ca", 2, 3), ("Fe", 2, 1), ("Mn", 4, 1)])
    mock_db = MagicMock(spec=StructureDB)
    mock_db.get_with_species.return_value = [parent]
    sp = StructurePredictor(cm, mock_db, table)
    predictions = list(sp.predict_structs(species, thresh=0.0, include_same=False))
    assert isinstance(predictions, list)


def test_predict_keyerror_lambda_table(predictor_coverage_data):
    db = predictor_coverage_data["db"]
    cm = predictor_coverage_data["cm"]
    table = predictor_coverage_data["table"]

    sp = StructurePredictor(cm, db, table)
    # Parent ca_fe_o3 has Fe4+ which is NOT in the lambda table.
    # cond_sub_probs("Ti4+").loc["Fe4+"] → KeyError → caught and continued
    species = [("Ca", 2), ("Ti", 4), ("O", -2)]
    predictions = list(sp.predict_structs(species, thresh=0.0, include_same=False))
    assert isinstance(predictions, list)


def test_predict_valueerror_non_neutral_mutation(predictor_coverage_data):
    db = predictor_coverage_data["db"]
    cm = predictor_coverage_data["cm"]
    table = predictor_coverage_data["table"]

    sp = StructurePredictor(cm, db, table)
    # Target K+, Cl-. DB has non-neutral Na+(x2) Cl-(x1).
    # Substituting Na+→K+ preserves non-neutrality → ValueError
    species = [("K", 1), ("Cl", -1)]
    predictions = list(sp.predict_structs(species, thresh=None, include_same=False))
    assert isinstance(predictions, list)


def test_nary_n1_parent_has_diff_species(predictor_coverage_data):
    db = predictor_coverage_data["db"]
    cm = predictor_coverage_data["cm"]
    table = predictor_coverage_data["table"]

    sp = StructurePredictor(cm, db, table)
    # NaCl in DB: querying for Na1+ returns NaCl which has Cl1- (the diff species)
    species = [("Na", 1), ("Cl", -1)]
    predictions = list(sp.nary_predict_structs(species, n_ary=1, include_same=False))
    assert isinstance(predictions, list)


def test_nary_n3_parent_has_all_diff_species(predictor_coverage_data):
    db = predictor_coverage_data["db"]
    cm = predictor_coverage_data["cm"]
    table = predictor_coverage_data["table"]

    sp = StructurePredictor(cm, db, table)
    # 4-species target. n_ary=3 queries for single species.
    # ba_ca_ti_o4 has all 4 → parent has all 3 diff species → continue
    species = [("Ba", 2), ("Ca", 2), ("Ti", 4), ("O", -2)]
    predictions = list(sp.nary_predict_structs(species, n_ary=3, include_same=False))
    assert isinstance(predictions, list)
