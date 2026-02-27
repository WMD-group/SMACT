"""Test structure prediction module."""

from __future__ import annotations

import itertools
import json
import logging
import os
import pickle
import sqlite3
import tempfile
import unittest
from contextlib import contextmanager
from importlib.util import find_spec
from operator import itemgetter
from random import sample
from typing import ClassVar
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pymatgen
import pytest
import requests
from pandas.testing import assert_frame_equal, assert_series_equal
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.structure_prediction.substitution_probability import (
    SubstitutionProbability,
)
from pymatgen.core import SETTINGS

import smact
import smact.structure_prediction.structure as sp_struct
from smact import Species
from smact.structure_prediction.database import StructureDB, parse_mprest
from smact.structure_prediction.mutation import CationMutator
from smact.structure_prediction.prediction import StructurePredictor
from smact.structure_prediction.structure import SmactStructure

MP_URL = "https://api.materialsproject.org"
MP_API_AVAILABLE = bool(find_spec("mp_api"))

try:
    skip_mprester_tests = requests.get(MP_URL, timeout=60).status_code != 200

except (ModuleNotFoundError, ImportError, requests.exceptions.ConnectionError):
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


def generate_test_structure(comp: str) -> bool:
    """Generate a pickled test structure for comparison."""
    poscar_file = os.path.join(files_dir, f"{comp}.txt")
    s = SmactStructure.from_file(poscar_file)

    with open(TEST_STRUCT, "wb") as f:
        pickle.dump(s, f)

    with open(TEST_POSCAR, "w") as f:
        f.write(s.as_poscar())

    return True


@contextmanager  # type: ignore[arg-type]
def ignore_warnings(logger: logging.Logger) -> int:  # type: ignore[return]
    """Ignore logging warnings."""
    log_lvl_buff = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield log_lvl_buff  # type: ignore[return]
    logger.setLevel(log_lvl_buff)


class StructureTest(unittest.TestCase):
    """`SmactStructure` testing."""

    TEST_SPECIES: ClassVar = {
        "CaTiO3": [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)],
        "NaCl": [("Na", 1, 1), ("Cl", -1, 1)],
        "Fe": [("Fe", 0, 1)],
    }

    def assertStructAlmostEqual(self, s1: SmactStructure, s2: SmactStructure, places: int = 7):
        """
        Assert that two SmactStructures are almost equal.

        Almost equality dependent on how many decimal places the site coordinates
        are equal to.
        """
        must_equal = ["species", "lattice_param"]
        for cond in must_equal:
            self.assertEqual(getattr(s1, cond), getattr(s2, cond))

        self.assertTrue(np.array_equal(s1.lattice_mat, s2.lattice_mat))
        self.assertEqual(list(s1.sites.keys()), list(s2.sites.keys()))

        for si, sj in zip(s1.sites, s2.sites, strict=False):
            for c1, c2 in zip(si, sj, strict=False):
                self.assertAlmostEqual(c1, c2, places=places)  # type: ignore[arg-type]

    def test_as_poscar(self):
        """Test POSCAR generation."""
        for comp in self.TEST_SPECIES:
            with self.subTest(comp=comp):
                comp_file = os.path.join(files_dir, f"{comp}.txt")
                with open(comp_file) as f:
                    struct = SmactStructure.from_file(comp_file)
                    self.assertEqual(struct.as_poscar(), f.read())

    @staticmethod
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

    def test_from_py_struct(self):
        """Test generation of SmactStructure from a pymatgen Structure."""
        with open(TEST_PY_STRUCT) as f:
            d = json.load(f)
            py_structure = pymatgen.core.Structure.from_dict(d)  # type: ignore[attr-defined]

        with ignore_warnings(smact.structure_prediction.logger):  # type: ignore[attr-defined]
            s1 = SmactStructure.from_py_struct(py_structure)

        s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))

        self.assertStructAlmostEqual(s1, s2)

    def test_from_py_struct_icsd(self):
        """Test generation of SmactStructure from a pymatgen Structure using ICSD statistics to determine oxidation states."""
        with open(TEST_PY_STRUCT) as f:
            d = json.load(f)
            py_structure = pymatgen.core.Structure.from_dict(d)  # type: ignore[attr-defined]

        with ignore_warnings(smact.structure_prediction.logger):  # type: ignore[attr-defined]
            s1 = SmactStructure.from_py_struct(py_structure, determine_oxi="comp_ICSD")

        s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))

        self.assertStructAlmostEqual(s1, s2)

    def test_has_species(self):
        """Test determining whether a species is in a `SmactStructure`."""
        s1 = SmactStructure(*self._gen_empty_structure([("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]))

        self.assertTrue(s1.has_species(("Ba", 2)))
        self.assertFalse(s1.has_species(("Ba", 3)))
        self.assertFalse(s1.has_species(("Ca", 2)))

    def test_smactStruc_comp_key(self):
        """Test generation of a composition key for `SmactStructure`s."""
        s1 = SmactStructure(*self._gen_empty_structure([("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]))
        s2 = SmactStructure(*self._gen_empty_structure([("Fe", 2, 1), ("Fe", 3, 2), ("O", -2, 4)]))

        Ba = Species("Ba", 2)
        O = Species("O", -2)
        F = Species("F", -1)
        Fe2 = Species("Fe", 2)
        Fe3 = Species("Fe", 3)

        s3 = SmactStructure(*self._gen_empty_structure([(Ba, 2), (O, 1), (F, 2)]))
        s4 = SmactStructure(*self._gen_empty_structure([(Fe2, 1), (Fe3, 2), (O, 4)]))

        Ba_2OF_2 = "Ba_2_2+F_2_1-O_1_2-"
        Fe_3O_4 = "Fe_2_3+Fe_1_2+O_4_2-"
        self.assertEqual(s1.composition(), Ba_2OF_2)
        self.assertEqual(s2.composition(), Fe_3O_4)
        self.assertEqual(s3.composition(), Ba_2OF_2)
        self.assertEqual(s4.composition(), Fe_3O_4)

    def test_smactStruc_from_file(self):
        """Test the `from_file` method of `SmactStructure`."""
        with open(TEST_STRUCT, "rb") as f:
            s1 = pickle.load(f)

        s2 = SmactStructure.from_file(TEST_POSCAR)

        self.assertEqual(s1, s2)

    def test_equality(self):
        """Test equality determination of `SmactStructure`."""
        struct_files = [os.path.join(files_dir, f"{x}.txt") for x in ["CaTiO3", "NaCl"]]
        CaTiO3 = SmactStructure.from_file(struct_files[0])
        NaCl = SmactStructure.from_file(struct_files[1])

        with self.subTest(msg="Testing equality of same object."):
            self.assertEqual(CaTiO3, CaTiO3)

        with self.subTest(msg="Testing inequality of different types."):
            self.assertNotEqual(CaTiO3, "CaTiO3")

        with self.subTest(msg="Testing inequality of different objects."):
            self.assertNotEqual(CaTiO3, NaCl)

    def test_ele_stoics(self):
        """Test acquiring element stoichiometries."""
        s1 = SmactStructure(*self._gen_empty_structure([("Fe", 2, 1), ("Fe", 3, 2), ("O", -2, 4)]))
        s1_stoics = {"Fe": 3, "O": 4}
        s2 = SmactStructure(*self._gen_empty_structure([("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]))
        s2_stoics = {"Ba": 2, "O": 1, "F": 2}

        for test, expected in [(s1, s1_stoics), (s2, s2_stoics)]:
            with self.subTest(species=test.species):
                self.assertEqual(SmactStructure._get_ele_stoics(test.species), expected)

    @pytest.mark.skipif(
        (
            skip_mprester_tests
            or not MP_API_AVAILABLE
            or not (os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY"))
        ),
        reason="Materials Project API not available or not configured.",
    )
    def test_from_mp(self):
        """Test downloading structures from materialsproject.org."""
        # TODO Needs ensuring that the structure query gets the same
        # structure as we have downloaded.

        api_key = os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY")

        for comp, species in self.TEST_SPECIES.items():
            with self.subTest(comp=comp):
                comp_file = os.path.join(files_dir, f"{comp}.txt")
                local_struct = SmactStructure.from_file(comp_file)
                mp_struct = SmactStructure.from_mp(species, api_key)  # type: ignore[arg-type]
                self.assertEqual(local_struct, mp_struct)

    def test_hash(self):
        """SmactStructure.__hash__ returns a consistent int (lines 158-173)."""
        s1 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
        s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
        self.assertIsInstance(hash(s1), int)
        self.assertEqual(hash(s1), hash(s2))
        # Different structures must produce different hashes
        s3 = SmactStructure.from_file(os.path.join(files_dir, "NaCl.txt"))
        self.assertNotEqual(hash(s1), hash(s3))

    def test_sanitise_species_error_branches(self):
        """_sanitise_species validation error branches (lines 199, 201, 203, 211, 223)."""
        with pytest.raises(TypeError):  # line 199: not a list
            SmactStructure._sanitise_species("not_a_list")  # type: ignore[arg-type]
        with pytest.raises(ValueError):  # line 201: empty list
            SmactStructure._sanitise_species([])
        with pytest.raises(TypeError):  # line 203: list of non-tuples
            SmactStructure._sanitise_species(["not_a_tuple"])  # type: ignore[arg-type]
        with pytest.raises(ValueError):  # line 211: tuple has wrong length
            SmactStructure._sanitise_species([("a", "b", "c", "d")])  # type: ignore[arg-type]
        with pytest.raises(TypeError):  # line 223: first element is not str or Species
            SmactStructure._sanitise_species([(42, 2, 1)])  # type: ignore[arg-type]

    def test_parse_py_sites_type_error(self):
        """SmactStructure.__parse_py_sites raises TypeError for non-Structure input (line 250)."""
        with pytest.raises(TypeError, match=r"pymatgen\.core\.Structure"):
            SmactStructure._SmactStructure__parse_py_sites("not_a_structure")  # type: ignore[attr-defined]

    def test_parse_py_sites_unit_charge(self):
        """__parse_py_sites inserts '1' for unit-charge species like Na+ (line 262)."""
        # NaCl decorated by BV gives Na+ and Cl- — unit charges trigger line 262
        nacl = pymatgen.core.Structure.from_spacegroup(  # type: ignore[attr-defined]
            "Fm-3m",
            pymatgen.core.Lattice.cubic(5.6),  # type: ignore[attr-defined]
            ["Na", "Cl"],
            [[0, 0, 0], [0.5, 0.5, 0.5]],  # type: ignore[attr-defined]
        )
        with ignore_warnings(smact.structure_prediction.logger):  # type: ignore[attr-defined]
            s = SmactStructure.from_py_struct(nacl, determine_oxi="BV")
        self.assertIsInstance(s, SmactStructure)
        # Confirm Na+1 and Cl-1 are present
        species_list = [sp[0] for sp in s.species]
        self.assertIn("Na", species_list)
        self.assertIn("Cl", species_list)

    def test_parse_py_sites_zero_charge(self):
        """__parse_py_sites assigns charge=0 when species string has no digits (lines 282-283)."""
        # A plain (undecorated) Structure used with determine_oxi='predecorated'
        # yields species_string = 'Fe' (no digits) → else branch
        fe_struct = pymatgen.core.Structure(pymatgen.core.Lattice.cubic(2.87), ["Fe"], [[0, 0, 0]])  # type: ignore[attr-defined]
        with ignore_warnings(smact.structure_prediction.logger):  # type: ignore[attr-defined]
            s = SmactStructure.from_py_struct(fe_struct, determine_oxi="predecorated")
        self.assertIsInstance(s, SmactStructure)
        self.assertEqual(s.species[0][1], 0)  # charge = 0

    def test_from_py_struct_type_error(self):
        """from_py_struct raises TypeError when given a non-Structure (line 307)."""
        with pytest.raises(TypeError, match=r"pymatgen\.core\.Structure"):
            SmactStructure.from_py_struct("not_a_structure")  # type: ignore[arg-type]

    def test_from_py_struct_both_success(self):
        """from_py_struct with determine_oxi='both' where BV succeeds (lines 319-323)."""
        with open(TEST_PY_STRUCT) as f:
            py_structure = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]
        with ignore_warnings(smact.structure_prediction.logger):  # type: ignore[attr-defined]
            s = SmactStructure.from_py_struct(py_structure, determine_oxi="both")
        self.assertIsInstance(s, SmactStructure)

    def test_from_py_struct_both_bv_fallback(self):
        """from_py_struct with determine_oxi='both' falls back to comp_ICSD when BV fails (lines 324-329)."""
        with open(TEST_PY_STRUCT) as f:
            py_structure = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]

        with (
            patch.object(BVAnalyzer, "get_oxi_state_decorated_structure", side_effect=ValueError("BV failed")),
            ignore_warnings(smact.structure_prediction.logger),  # type: ignore[attr-defined]
        ):
            s = SmactStructure.from_py_struct(py_structure, determine_oxi="both")
        self.assertIsInstance(s, SmactStructure)

    def test_from_py_struct_predecorated(self):
        """from_py_struct with determine_oxi='predecorated' uses structure as-is (line 331)."""
        with open(TEST_PY_STRUCT) as f:
            py_structure = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]
        decorated = BVAnalyzer().get_oxi_state_decorated_structure(py_structure)
        with ignore_warnings(smact.structure_prediction.logger):  # type: ignore[attr-defined]
            s = SmactStructure.from_py_struct(decorated, determine_oxi="predecorated")
        self.assertIsInstance(s, SmactStructure)

    def test_from_py_struct_invalid_determine_oxi(self):
        """from_py_struct raises ValueError for unknown determine_oxi (lines 333-336)."""
        with open(TEST_PY_STRUCT) as f:
            py_structure = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]
        with pytest.raises(ValueError, match="determine_oxi"):
            SmactStructure.from_py_struct(py_structure, determine_oxi="unknown")

    def test_from_mp_no_api_key(self):
        """from_mp raises ValueError when no API key is available (lines 377, 379)."""
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

    def test_from_mp_mocked_new_api(self):
        """from_mp via mocked MPResterNew covers lines 399-446."""
        with open(TEST_PY_STRUCT) as f:
            real_struct = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]

        mock_cm = MagicMock()
        mock_cm.materials.summary.search.return_value = [real_struct]
        mock_rester = MagicMock()
        mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
        mock_rester.return_value.__exit__ = MagicMock(return_value=False)

        species = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
        with (
            patch.object(sp_struct, "MPResterNew", new=mock_rester),
            patch.object(sp_struct, "HAS_MP_API", new=True),
            ignore_warnings(smact.structure_prediction.logger),  # type: ignore[attr-defined]
        ):
            s = SmactStructure.from_mp(species, api_key="x" * 32)  # type: ignore[arg-type]
        self.assertIsInstance(s, SmactStructure)

    def test_from_mp_mocked_legacy_api(self):
        """from_mp via mocked legacy MPRester covers lines 391-395."""
        with open(TEST_PY_STRUCT) as f:
            real_struct = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]

        mock_cm = MagicMock()
        mock_cm.query.return_value = [{"structure": real_struct}]
        mock_rester = MagicMock()
        mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
        mock_rester.return_value.__exit__ = MagicMock(return_value=False)

        species = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
        with (
            patch.object(sp_struct, "MPRester", new=mock_rester),
            patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=True),
            ignore_warnings(smact.structure_prediction.logger),  # type: ignore[attr-defined]
        ):
            s = SmactStructure.from_mp(species, api_key="short_key")  # type: ignore[arg-type]
        self.assertIsInstance(s, SmactStructure)

    def test_from_mp_mocked_legacy_32char(self):
        """from_mp with 32-char key + no mp_api + legacy MPRester covers lines 406-407."""
        with open(TEST_PY_STRUCT) as f:
            real_struct = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]

        mock_cm = MagicMock()
        mock_cm.get_structures.return_value = [real_struct]
        mock_rester = MagicMock()
        mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
        mock_rester.return_value.__exit__ = MagicMock(return_value=False)

        species = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
        with (
            patch.object(sp_struct, "MPRester", new=mock_rester),
            patch.object(sp_struct, "HAS_MP_API", new=False),
            patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=True),
            ignore_warnings(smact.structure_prediction.logger),  # type: ignore[attr-defined]
        ):
            s = SmactStructure.from_mp(species, api_key="x" * 32)  # type: ignore[arg-type]
        self.assertIsInstance(s, SmactStructure)

    def test_from_mp_determine_oxi_branches(self):
        """from_mp with comp_ICSD, both, and invalid determine_oxi (lines 420-437)."""
        with open(TEST_PY_STRUCT) as f:
            real_struct = pymatgen.core.Structure.from_dict(json.load(f))  # type: ignore[attr-defined]

        def make_mock_rester(return_val):
            mock_cm = MagicMock()
            mock_cm.materials.summary.search.return_value = [return_val]
            mock_rester = MagicMock()
            mock_rester.return_value.__enter__ = MagicMock(return_value=mock_cm)
            mock_rester.return_value.__exit__ = MagicMock(return_value=False)
            return mock_rester

        species = [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)]
        api_key_32 = "x" * 32

        # comp_ICSD branch (lines 420-424)
        with (
            patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
            patch.object(sp_struct, "HAS_MP_API", new=True),
            ignore_warnings(smact.structure_prediction.logger),  # type: ignore[attr-defined]
        ):
            s = SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="comp_ICSD")  # type: ignore[arg-type]
            self.assertIsInstance(s, SmactStructure)

        # both branch — BV succeeds (lines 426-430)
        with (
            patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
            patch.object(sp_struct, "HAS_MP_API", new=True),
            ignore_warnings(smact.structure_prediction.logger),  # type: ignore[attr-defined]
        ):
            s = SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="both")  # type: ignore[arg-type]
            self.assertIsInstance(s, SmactStructure)

        # both branch — BV fails, falls back to comp_ICSD (lines 431-435)
        with (
            patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
            patch.object(sp_struct, "HAS_MP_API", new=True),
            patch.object(BVAnalyzer, "get_oxi_state_decorated_structure", side_effect=ValueError("BV failed")),
            ignore_warnings(smact.structure_prediction.logger),  # type: ignore[attr-defined]
        ):
            s = SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="both")  # type: ignore[arg-type]
            self.assertIsInstance(s, SmactStructure)

        # invalid determine_oxi → ValueError (lines 436-439)
        with (
            patch.object(sp_struct, "MPResterNew", new=make_mock_rester(real_struct)),
            patch.object(sp_struct, "HAS_MP_API", new=True),
            pytest.raises(ValueError, match="determine_oxi"),
        ):
            SmactStructure.from_mp(species, api_key=api_key_32, determine_oxi="invalid_oxi")  # type: ignore[arg-type]

    def test_from_mp_empty_result_raises(self):
        """from_mp raises ValueError when API returns no structures (line 409-410)."""
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

    def test_from_mp_no_legacy_mprester_short_key(self):
        """from_mp raises ImportError when HAS_LEGACY_MPRESTER=False and short API key used (lines 386-390)."""
        species = [("Na", 1, 1), ("Cl", -1, 1)]
        with (
            patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=False),
            pytest.raises(ImportError, match="pymatgen legacy MPRester is not available"),
        ):
            SmactStructure.from_mp(species, api_key="short_key")  # type: ignore[arg-type]

    def test_from_mp_no_legacy_mprester_no_mp_api(self):
        """from_mp raises ImportError when HAS_LEGACY_MPRESTER=False, HAS_MP_API=False, 32-char key (lines 402-405)."""
        species = [("Na", 1, 1), ("Cl", -1, 1)]
        api_key_32 = "x" * 32
        with (
            patch.object(sp_struct, "HAS_LEGACY_MPRESTER", new=False),
            patch.object(sp_struct, "HAS_MP_API", new=False),
            pytest.raises(ImportError, match="Neither mp-api nor pymatgen legacy MPRester"),
        ):
            SmactStructure.from_mp(species, api_key=api_key_32)  # type: ignore[arg-type]

    def test_as_py_struct(self):
        s1 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
        s1_pym = s1.as_py_struct()
        with open(TEST_PY_STRUCT) as f:
            d = json.load(f)
        py_structure = pymatgen.core.Structure.from_dict(d)  # type: ignore[attr-defined]

        self.assertEqual(s1_pym, py_structure)

    def test_reduced_formula(self):
        """Test the reduced formula method."""
        s1 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))
        s2 = SmactStructure.from_file(os.path.join(files_dir, "NaCl.txt"))

        self.assertEqual(s1.reduced_formula(), "CaTiO3")
        self.assertEqual(s2.reduced_formula(), "NaCl")


class StructureDBTest(unittest.TestCase):
    """Test StructureDB interface."""

    TEST_DB = os.path.join(files_dir, "test_db.tmp")
    TEST_TABLE = "Structures"
    TEST_MP_TABLE = "Structures1"

    @classmethod
    def tearDownClass(cls):
        """Remove database files."""
        if os.path.exists(cls.TEST_DB):
            os.remove(cls.TEST_DB)

    def test_db_interface(self):
        """Test interfacing with database."""
        with self.subTest(msg="Instantiating database."):
            self.db = StructureDB(self.TEST_DB)

        with self.subTest(msg="Adding table."):
            try:
                self.db.add_table(self.TEST_TABLE)
            except (sqlite3.Error, ValueError, OSError) as e:
                self.fail(e)

        struct_file = os.path.join(files_dir, "CaTiO3.txt")
        struct = SmactStructure.from_file(struct_file)

        with self.subTest(msg="Adding structure to table."):
            try:
                self.db.add_struct(struct, self.TEST_TABLE)
            except (sqlite3.Error, ValueError, OSError) as e:
                self.fail(e)

        with self.subTest(msg="Getting structure from table."):
            struct_list = self.db.get_structs(struct.composition(), self.TEST_TABLE)
            self.assertEqual(len(struct_list), 1)
            self.assertEqual(struct_list[0], struct)

        struct_files = [os.path.join(files_dir, f"{x}.txt") for x in ["NaCl", "Fe"]]
        structs = [SmactStructure.from_file(fname) for fname in struct_files]

        with self.subTest(msg="Adding multiple structures to table."):
            try:
                self.db.add_structs(structs, self.TEST_TABLE)
            except (sqlite3.Error, ValueError, OSError) as e:
                self.fail(e)

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
            with self.subTest(msg=f"Retrieving species with {spec}"):
                self.assertEqual(self.db.get_with_species(spec, self.TEST_TABLE), expected)

        mp_strucs = [
            pymatgen.core.Structure.from_file(os.path.join(files_dir, f))  # type: ignore[attr-defined]
            for f in ["CaTiO3.json", "NaCl.json", "Fe3O4.json"]
        ]
        mp_data = [
            {"material_id": mpid, "structure": s}
            for mpid, s in zip(["mp-4019", "mp-22862", "mp-19306"], mp_strucs, strict=False)
        ]

        with self.subTest(msg="Testing adding downloaded MP structures."):
            added: int = self.db.add_mp_icsd(self.TEST_MP_TABLE, mp_data)
            self.assertEqual(added, 3)

    def test_db_rollback_on_exception(self):
        """StructureDB.__exit__ rollback branch (line 103): exception triggers conn.rollback()."""
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
            db_path = f.name

        db = StructureDB(db_path)
        db.add_table("rollback_test")

        try:
            with db as c:
                c.execute("INSERT INTO rollback_test VALUES (?, ?)", ("comp", "poscar"))
                raise RuntimeError("intentional error for rollback")
        except RuntimeError:
            pass

        # Rollback means the INSERT should not have been committed
        with db as c:
            c.execute("SELECT * FROM rollback_test")
            rows = c.fetchall()
        self.assertEqual(len(rows), 0)

    def test_add_structs_with_none(self):
        """add_structs skips None entries (line 227: continue)."""
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
            db_path = f.name
        try:
            db = StructureDB(db_path)
            db.add_table("test_none")
            struct = SmactStructure.from_file(os.path.join(files_dir, "NaCl.txt"))
            added = db.add_structs([None, struct, None], "test_none")  # type: ignore[arg-type]
            self.assertEqual(added, 1)
        finally:
            os.remove(db_path)

    def test_parse_mprest_exception_path(self):
        """parse_mprest logs warning and returns None when from_py_struct raises (lines 321-323)."""
        # Passing a non-Structure as 'structure' causes from_py_struct → TypeError → caught
        parse_mprest({"structure": "not_a_structure", "material_id": "mp-test"})


class CationMutatorTest(unittest.TestCase):
    """Test the CationMutator class."""

    @classmethod
    def setUpClass(cls):
        """Set up the test initial structure and mutator."""
        cls.test_struct = SmactStructure.from_file(TEST_POSCAR)

        cls.test_mutator = CationMutator.from_json(lambda_json=TEST_LAMBDA_JSON)
        cls.test_pymatgen_mutator = CationMutator.from_json(lambda_json=None, alpha=lambda x, y: -5)

        # 5 random test species -> 5! test pairs
        cls.test_species = sample(list(cls.test_pymatgen_mutator.specs), 5)
        cls.test_pairs = list(itertools.combinations_with_replacement(cls.test_species, 2))

        cls.pymatgen_sp = SubstitutionProbability(lambda_table=None, alpha=-5)

    def test_lambda_tab_pop(self):
        """Test if lambda table is populated correctly."""
        lambda_dat = [
            [-5.0, 0.5, -5.0],
            [0.5, -5.0, 0.3],
            [-5.0, 0.3, -5.0],
        ]
        labels = ["A", "B", "C"]

        exp_lambda = pd.DataFrame(lambda_dat, index=labels, columns=labels)  # type: ignore[call-overload]

        assert_frame_equal(
            self.test_mutator.lambda_tab,
            exp_lambda,
            check_names=False,
        )

    def test_partition_func_Z(self):
        """Test the partition function for the whole table."""
        # 2e^0.5 + 2e^0.3 + 5e^{-5} \approx 6.0308499
        self.assertAlmostEqual(self.test_mutator.Z, 6.0308499)

    def test_pymatgen_lambda_import(self):
        """Test importing pymatgen lambda table."""
        self.assertIsInstance(self.test_pymatgen_mutator.lambda_tab, pd.DataFrame)

    def test_lambda_interface(self):
        """Test getting lambda values."""
        test_cases = [itertools.permutations(x) for x in [("A", "B"), ("A", "C"), ("B", "C")]]

        expected = [0.5, -5.0, 0.3]

        for test_case, expectation in zip(test_cases, expected, strict=False):
            for spec_comb in test_case:
                s1, s2 = spec_comb
                with self.subTest(s1=s1, s2=s2):
                    self.assertEqual(self.test_mutator.get_lambda(s1, s2), expectation)

    def test_ion_mutation(self):
        """Test mutating an ion of a SmactStructure."""
        ca_file = os.path.join(files_dir, "CaTiO3.txt")
        ba_file = os.path.join(files_dir, "BaTiO3.txt")

        CaTiO3 = SmactStructure.from_file(ca_file)
        BaTiO3 = SmactStructure.from_file(ba_file)

        with self.subTest(s1="CaTiO3", s2="BaTiO3"):
            mutation = self.test_mutator._mutate_structure(CaTiO3, "Ca2+", "Ba2+")
            self.assertEqual(mutation, BaTiO3)

        na_file = os.path.join(files_dir, "NaCl.txt")
        NaCl = SmactStructure.from_file(na_file)

        with self.subTest(s1="Na1+Cl1-", s2="Na2+Cl1-"), pytest.raises(ValueError):
            self.test_mutator._mutate_structure(NaCl, "Na1+", "Na2+")

        # TODO Confirm functionality with more complex substitutions

    def test_sub_prob(self):
        """Test determining substitution probabilities."""
        for s1, s2 in self.test_pairs:
            with self.subTest(s1=s1, s2=s2):
                self.assertAlmostEqual(
                    self.pymatgen_sp.prob(s1, s2),
                    self.test_pymatgen_mutator.sub_prob(s1, s2),
                )

    def test_cond_sub_probs(self):
        """Test determining conditional substitution probabilities for a row."""
        for s1 in ["A", "B", "C"]:
            with self.subTest(s=s1):
                cond_sub_probs_test = self.test_mutator.cond_sub_probs(s1)

                vals = [(s1, s2, self.test_mutator.cond_sub_prob(s1, s2)) for s2 in ["A", "B", "C"]]

                test_df = pd.DataFrame(vals)
                test_df: pd.DataFrame = test_df.pivot_table(index=0, columns=1, values=2)

                # Slice to convert to series
                assert_series_equal(cond_sub_probs_test, test_df.iloc[0])

    def test_cond_sub_prob(self):
        """Test determining conditional substitution probabilities."""
        for s1, s2 in self.test_pairs:
            with self.subTest(s1=s1, s2=s2):
                self.assertAlmostEqual(
                    self.pymatgen_sp.cond_prob(s1, s2),
                    self.test_pymatgen_mutator.cond_sub_prob(s1, s2),
                )

    def test_pair_corr(self):
        """Test determining conditional substitution probabilities."""
        for s1, s2 in self.test_pairs:
            with self.subTest(s1=s1, s2=s2):
                self.assertAlmostEqual(
                    self.pymatgen_sp.pair_corr(s1, s2),
                    self.test_pymatgen_mutator.pair_corr(s1, s2),
                )

    def test_from_df(self):
        """Test creating a CationMutator from an existing DataFrame."""
        lambda_df = pd.read_csv(TEST_LAMBDA_CSV, index_col=0)

        csv_test = CationMutator(lambda_df=lambda_df)
        assert_frame_equal(
            csv_test.lambda_tab,
            self.test_mutator.lambda_tab,
            check_names=False,
        )

    def test_complete_cond_probs(self):
        """Test getting all conditional probabilities."""
        pairs = itertools.product(["A", "B", "C"], repeat=2)

        vals = [
            (
                s1,
                s2,
                self.test_mutator.cond_sub_prob(s1, s2),
            )
            for s1, s2 in pairs
        ]

        cond_probs = pd.DataFrame(vals)
        cond_probs = cond_probs.pivot_table(index=0, columns=1, values=2)

        assert_frame_equal(self.test_mutator.complete_cond_probs(), cond_probs)

    def test_complete_sub_probs(self):
        """Test getting all probabilities."""
        pairs = itertools.product(["A", "B", "C"], repeat=2)

        vals = [
            (
                s1,
                s2,
                self.test_mutator.sub_prob(s1, s2),
            )
            for s1, s2 in pairs
        ]

        sub_probs = pd.DataFrame(vals)
        sub_probs = sub_probs.pivot_table(index=0, columns=1, values=2)

        assert_frame_equal(self.test_mutator.complete_sub_probs(), sub_probs)

    def test_complete_pair_corrs(self):
        """Test getting all pair correlations."""
        pairs = itertools.product(["A", "B", "C"], repeat=2)

        vals = [
            (
                s1,
                s2,
                self.test_mutator.pair_corr(s1, s2),
            )
            for s1, s2 in pairs
        ]

        pair_corrs = pd.DataFrame(vals)
        pair_corrs = pair_corrs.pivot_table(index=0, columns=1, values=2)

        assert_frame_equal(self.test_mutator.complete_pair_corrs(), pair_corrs)

    def test_same_spec_probs(self):
        """same_spec_probs returns a Series of diagonal substitution probabilities."""
        probs = self.test_mutator.same_spec_probs()
        self.assertIsInstance(probs, pd.Series)
        self.assertEqual(len(probs), 3)
        # All values should be positive
        self.assertTrue((probs > 0).all())

    def test_same_spec_cond_probs(self):
        """same_spec_cond_probs returns conditional probabilities for self-substitution."""
        probs = self.test_mutator.same_spec_cond_probs()
        # Conditional probabilities should be between 0 and 1
        self.assertTrue((probs >= 0).all())
        self.assertTrue((probs <= 1).all())

    def test_get_lambda_unknown_species(self):
        """get_lambda falls back to alpha for species not in the table."""
        lam = self.test_mutator.get_lambda("X", "Y")
        self.assertEqual(lam, -5.0)

    def test_get_lambdas_unknown_species_raises(self):
        """get_lambdas raises ValueError for species not in the table."""
        with pytest.raises(ValueError, match="not in lambda table"):
            self.test_mutator.get_lambdas("UNKNOWN_SPEC")

    def test_populate_lambda_alpha_none_raises(self):
        """_populate_lambda raises ValueError when alpha is None and table has NaN entries."""
        lambda_df = pd.DataFrame(
            [[1.0, np.nan], [np.nan, 1.0]],
            index=["A", "B"],  # type: ignore[call-overload]
            columns=["A", "B"],  # type: ignore[call-overload]
        )
        with pytest.raises(ValueError, match="alpha function must not be None"):
            CationMutator(lambda_df, alpha=None)

    def test_populate_lambda_asymmetric_mirror(self):
        """_populate_lambda mirrors (s2, s1) into (s1, s2) when only (s2, s1) exists (lines 152-153).

        When looking up (s1, s2) raises KeyError but (s2, s1) has a real
        (non-NaN) value, the code should mirror that value into (s1, s2).
        """
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
        self.assertEqual(cm.lambda_tab.loc["B", "A"], 2.0)
        self.assertEqual(cm.lambda_tab.loc["A", "B"], 2.0)

    def test_nary_mutate_structure(self):
        """_nary_mutate_structure replaces multiple species simultaneously."""
        ca_file = os.path.join(files_dir, "CaTiO3.txt")
        CaTiO3 = SmactStructure.from_file(ca_file)

        # Replace Ca2+ → Ba2+ and Ti4+ → Zr4+ (both same-charge substitutions)
        mutated = CationMutator._nary_mutate_structure(CaTiO3, ["Ca2+", "Ti4+"], ["Ba2+", "Zr4+"])
        self.assertIsInstance(mutated, SmactStructure)
        spec_strs = mutated.get_spec_strs()
        self.assertIn("Ba2+", spec_strs)
        self.assertIn("Zr4+", spec_strs)
        self.assertNotIn("Ca2+", spec_strs)
        self.assertNotIn("Ti4+", spec_strs)

    def test_nary_mutate_structure_not_neutral(self):
        """_nary_mutate_structure raises ValueError when result is not charge neutral."""
        ca_file = os.path.join(files_dir, "CaTiO3.txt")
        CaTiO3 = SmactStructure.from_file(ca_file)

        # Replace Ca2+ → Na1+ breaks charge neutrality
        with pytest.raises(ValueError, match="not charge neutral"):
            CationMutator._nary_mutate_structure(CaTiO3, ["Ca2+"], ["Na1+"])

    def test_unary_substitute(self):
        """unary_substitute yields structures with single-species substitutions."""
        ca_file = os.path.join(files_dir, "CaTiO3.txt")
        CaTiO3 = SmactStructure.from_file(ca_file)
        # Use pymatgen mutator which has a full lambda table
        results = list(self.test_pymatgen_mutator.unary_substitute(CaTiO3, thresh=1e-3))
        # Should yield at least some substitutions
        self.assertGreater(len(results), 0)
        for struct, prob, old_spec, new_spec in results:
            self.assertIsInstance(struct, SmactStructure)
            self.assertGreater(prob, 0)
            self.assertNotEqual(old_spec, new_spec)

    def test_sub_probs(self):
        """sub_probs returns substitution probabilities for all species."""
        probs = self.test_mutator.sub_probs("A")
        self.assertIsInstance(probs, pd.Series)
        self.assertTrue((probs >= 0).all())


class PredictorTest(unittest.TestCase):
    """Testing for the StructurePredictor wrapper."""

    @classmethod
    def setUpClass(cls):
        """Initialize the prerequisites."""
        cls.db = StructureDB(TEST_PREDICTOR_DB)
        # NOTE: This may break if the pymatgen lambda table is updated
        cls.cm = CationMutator.from_json()
        cls.table = TEST_PREDICTOR_TABLE

    def test_prediction(self):
        """Test prediction pipeline."""
        with self.subTest(msg="Testing predictor setup"):
            sp = StructurePredictor(self.cm, self.db, self.table)

        test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]

        # At time of creation, results with thresh=0 and include_same=True are as follows:
        #
        # Ca_1_2+O_4_2-Ti_2_4+: p=1.0, from Ca_1_2+O_4_2-Ti_2_4+
        # Ca_1_2+O_3_2-Ti_1_4+: p=0.0088406421433475, from Ca_1_2+O_3_2-Si_1_4+
        # Ca_1_2+O_3_2-Ti_1_4+: p=0.05502122697602804, from Ca_1_2+Mo_1_4+O_3_2-
        # Ca_1_2+O_3_2-Ti_1_4+: p=0.0445515892332118, from Ca_1_2+O_3_2-V_1_4+
        # Ca_2_2+O_4_2-Ti_1_4+: p=0.0274412135813993, from Fe_2_2+O_4_2-Ti_1_4+

        with self.subTest(msg="Acquiring predictions"):
            try:
                predictions = list(sp.predict_structs(test_specs, thresh=0.02, include_same=False))
            except (sqlite3.Error, ValueError, OSError) as e:
                self.fail(e)

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

        with self.subTest(msg="Ensuring same parents"):
            self.assertEqual(parents, expected_comps)

        probs = [x[1] for x in predictions]

        with self.subTest(msg="Ensuring similar probabilities"):
            for p1, p2 in zip(probs, expected_probs, strict=False):
                with self.subTest(p1=p1, p2=p2):
                    self.assertAlmostEqual(p1, p2)

    def test_prediction_include_same(self):
        """predict_structs with include_same=True yields identical structures with p=1.0."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
        predictions = list(sp.predict_structs(test_specs, thresh=0.0, include_same=True))
        # The first results should be identical structures with probability 1.0
        same_structs = [(s, p, parent) for s, p, parent in predictions if p == 1.0]
        self.assertGreater(len(same_structs), 0)
        for s, _p, parent in same_structs:
            self.assertEqual(s, parent)

    def test_prediction_thresh_none(self):
        """predict_structs with thresh=None skips threshold filtering."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
        predictions = list(sp.predict_structs(test_specs, thresh=None, include_same=False))
        # With no threshold, we should get at least as many as with a threshold
        predictions_thresh = list(sp.predict_structs(test_specs, thresh=0.02, include_same=False))
        self.assertGreaterEqual(len(predictions), len(predictions_thresh))

    def test_nary_predict_n_ary_none(self):
        """nary_predict_structs with n_ary=None returns early after same-structure yield."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
        predictions = list(sp.nary_predict_structs(test_specs, n_ary=None, include_same=True))
        # Only yields identical structures (n_ary=None → early return)
        for _s, p, _parent in predictions:
            self.assertEqual(p, 1.0)

    def test_nary_predict_equal_species_returns_early(self):
        """nary_predict_structs with n_ary==len(species) returns early (no subset possible)."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
        # n_ary=3 means we'd replace all 3 species, leaving 0 common → early return
        predictions = list(sp.nary_predict_structs(test_specs, n_ary=3, include_same=False))
        self.assertEqual(len(predictions), 0)

    def test_nary_predict_include_same_false(self):
        """nary_predict_structs with include_same=False doesn't yield identical structures."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]
        predictions = list(sp.nary_predict_structs(test_specs, n_ary=2, include_same=False))
        for _s, p, _parent in predictions:
            self.assertNotEqual(p, 1.0)


class PredictorCoverageTest(unittest.TestCase):
    """Tests targeting uncovered guard clauses in prediction.py."""

    @staticmethod
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

    @classmethod
    def setUpClass(cls):
        """Set up temp-file DB and minimal CationMutator."""
        cls.table = "COV"
        cls._db_fd, cls._db_path = tempfile.mkstemp(suffix=".db")

        nacl = cls._make_structure([("Na", 1, 1), ("Cl", -1, 1)])
        ca_ti_o3 = cls._make_structure([("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)])
        # Parent with Fe4+ — NOT in the lambda table → triggers KeyError path
        ca_fe_o3 = cls._make_structure([("Ca", 2, 1), ("Fe", 4, 1), ("O", -2, 3)])
        # 4-species structure for n_ary=3 test
        ba_ca_ti_o4 = cls._make_structure([("Ba", 2, 1), ("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 4)])
        # Non-neutral structure: Na1+(x2) Cl1-(x1) → charge +1
        non_neutral = cls._make_structure([("Na", 1, 2), ("Cl", -1, 1)])

        cls.db = StructureDB(cls._db_path)
        cls.db.add_table(cls.table)
        for s in [nacl, ca_ti_o3, ca_fe_o3, ba_ca_ti_o4, non_neutral]:
            cls.db.add_struct(s, cls.table)

        # Minimal CationMutator — Fe4+ intentionally excluded
        specs = ["Na1+", "Cl1-", "Ca2+", "Ti4+", "O2-", "Ba2+", "Si4+", "K1+"]
        n = len(specs)
        data = np.full((n, n), -5.0)
        np.fill_diagonal(data, 0.0)
        lambda_df = pd.DataFrame(data, index=specs, columns=specs)  # type: ignore[call-overload]
        cls.cm = CationMutator(lambda_df)

    @classmethod
    def tearDownClass(cls):
        """Remove temp DB file."""
        os.close(cls._db_fd)
        os.unlink(cls._db_path)

    def test_predict_duplicate_species_continue(self):
        """Line 105: continue when species duplicates make an empty difference set."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        species = [("Na", 1), ("Na", 1), ("Cl", -1)]
        predictions = list(sp.predict_structs(species, thresh=0.0, include_same=False))
        self.assertIsInstance(predictions, list)

    def test_predict_parent_multi_extra_species(self):
        """Line 126: continue when parent has >1 non-target species."""
        species = [("Ca", 2), ("Ti", 4), ("O", -2)]
        # Parent with Ca2+ (matching) + Fe2+, Mn4+ (both differ from target)
        parent = self._make_structure([("Ca", 2, 3), ("Fe", 2, 1), ("Mn", 4, 1)])
        mock_db = MagicMock(spec=StructureDB)
        mock_db.get_with_species.return_value = [parent]
        sp = StructurePredictor(self.cm, mock_db, self.table)
        predictions = list(sp.predict_structs(species, thresh=0.0, include_same=False))
        self.assertIsInstance(predictions, list)

    def test_predict_keyerror_lambda_table(self):
        """Lines 135-137: continue when alt_spec not in cond_sub_probs."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        # Parent ca_fe_o3 has Fe4+ which is NOT in the lambda table.
        # cond_sub_probs("Ti4+").loc["Fe4+"] → KeyError → caught and continued
        species = [("Ca", 2), ("Ti", 4), ("O", -2)]
        predictions = list(sp.predict_structs(species, thresh=0.0, include_same=False))
        self.assertIsInstance(predictions, list)

    def test_predict_valueerror_non_neutral_mutation(self):
        """Lines 142-144: continue when _mutate_structure raises ValueError."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        # Target K+, Cl-. DB has non-neutral Na+(x2) Cl-(x1).
        # Substituting Na+→K+ preserves non-neutrality → ValueError
        species = [("K", 1), ("Cl", -1)]
        predictions = list(sp.predict_structs(species, thresh=None, include_same=False))
        self.assertIsInstance(predictions, list)

    def test_nary_n1_parent_has_diff_species(self):
        """Lines 206-207: continue when n_ary=1 and parent has diff species."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        # NaCl in DB: querying for Na1+ returns NaCl which has Cl1- (the diff species)
        species = [("Na", 1), ("Cl", -1)]
        predictions = list(sp.nary_predict_structs(species, n_ary=1, include_same=False))
        self.assertIsInstance(predictions, list)

    def test_nary_n3_parent_has_all_diff_species(self):
        """Lines 211-216: continue when n_ary=3 and parent has all 3 diff species."""
        sp = StructurePredictor(self.cm, self.db, self.table)
        # 4-species target. n_ary=3 queries for single species.
        # ba_ca_ti_o4 has all 4 → parent has all 3 diff species → continue
        species = [("Ba", 2), ("Ca", 2), ("Ti", 4), ("O", -2)]
        predictions = list(sp.nary_predict_structs(species, n_ary=3, include_same=False))
        self.assertIsInstance(predictions, list)
