from __future__ import annotations

import os
import shutil
import sys
import unittest
from importlib.util import find_spec

import pandas as pd
import pytest
import requests
from pymatgen.core import SETTINGS, Composition

from smact import Element
from smact.screening import smact_filter
from smact.utils.composition import comp_maker, formula_maker, parse_formula
from smact.utils.crystal_space import generate_composition_with_smact
from smact.utils.oxidation import ICSD24OxStatesFilter

MP_URL = "https://api.materialsproject.org"
MP_API_AVAILABLE = bool(find_spec("mp_api"))

try:
    skip_mprester_tests = requests.get(MP_URL, timeout=60).status_code != 200

except (ModuleNotFoundError, ImportError, requests.exceptions.ConnectionError):
    # Skip all MPRester tests if some downstream problem on the website, mp-api or whatever.
    skip_mprester_tests = True


class TestComposition(unittest.TestCase):
    """Test composition utilities"""

    def setUp(self) -> None:
        self.mock_filter_output = [
            (("Fe", "O"), (2, -2), (1, 1)),
            (("Fe", "O"), (1, 1)),
            (("Fe", "Fe", "O"), (2, 3, -2), (1, 2, 4)),
        ]
        self.smact_filter_output = smact_filter(
            els=[Element("Li"), Element("Ge"), Element("P"), Element("S")],
            stoichs=[[10], [1], [2], [12]],
        )

    def test_parse_formula(self):
        """Test the parse_formula function"""

        formulas = ["Li10GeP2S12", "Mg0.5O0.5", "CaMg(CO3)2"]

        LGPS = parse_formula(formulas[0])
        self.assertIsInstance(LGPS, dict)
        for el_sym, ammt in LGPS.items():
            self.assertIsInstance(el_sym, str)
            self.assertIsInstance(ammt, float)
        self.assertEqual(LGPS["Li"], 10)
        self.assertEqual(LGPS["Ge"], 1)
        self.assertEqual(LGPS["P"], 2)
        self.assertEqual(LGPS["S"], 12)

        MgO = parse_formula(formulas[1])
        self.assertIsInstance(MgO, dict)
        self.assertEqual(MgO["Mg"], 0.5)
        self.assertEqual(MgO["O"], 0.5)

        dolomite = parse_formula(formulas[2])
        self.assertIsInstance(dolomite, dict)
        self.assertEqual(dolomite["Ca"], 1)
        self.assertEqual(dolomite["Mg"], 1)
        self.assertEqual(dolomite["C"], 2)
        self.assertEqual(dolomite["O"], 6)

    def test_comp_maker(self):
        """Test the comp_maker function"""
        comp1 = comp_maker(self.mock_filter_output[0])
        comp2 = comp_maker(self.mock_filter_output[1])
        comp3 = comp_maker(self.mock_filter_output[2])
        comp4 = comp_maker(self.smact_filter_output[1])
        for comp in [comp1, comp2, comp3, comp4]:
            self.assertIsInstance(comp, Composition)
        self.assertEqual(Composition("FeO"), comp2)
        self.assertEqual(Composition({"Fe2+": 1, "O2-": 1}), comp1)
        self.assertEqual(Composition({"Fe2+": 1, "Fe3+": 2, "O2-": 4}), comp3)
        self.assertEqual(Composition({"Li+": 10, "Ge4+": 1, "P5+": 2, "S2-": 12}), comp4)

    def test_formula_maker(self):
        """Test the formula_maker function"""
        form1 = formula_maker(self.mock_filter_output[0])
        form2 = formula_maker(self.mock_filter_output[1])
        form3 = formula_maker(self.mock_filter_output[2])
        form4 = formula_maker(self.smact_filter_output[1])
        self.assertEqual(form1, "FeO")
        self.assertEqual(form2, "FeO")
        self.assertEqual(form1, form2)
        self.assertEqual(form3, "Fe3O4")
        self.assertEqual(form4, "Li10Ge(PS6)2")


class TestCrystalSpace(unittest.TestCase):
    """Test utility functions associated with Crystal Space."""

    def test_convert_formula(self):
        combinations = [Element("Li"), Element("O")]
        expected_formulas = ["Li2O2", "LiO2", "Li2O", "Li2O2"]
        compounds = generate_composition_with_smact.convert_formula(
            combinations=combinations, num_elements=2, max_stoich=2
        )
        self.assertListEqual(expected_formulas, compounds)

    def test_generate_composition_with_smact(self):
        save_dir = "data/binary/df_binary_label.pkl"
        oxidation_states_sets = ["smact14", "icsd24"]
        oxidation_states_sets_dict = {
            "smact14": {"smact_allowed": 388},
            "icsd24": {"smact_allowed": 342},
        }
        for ox_states in oxidation_states_sets:
            with self.subTest(ox_states=ox_states):
                smact_df = generate_composition_with_smact.generate_composition_with_smact(
                    num_elements=2,
                    max_stoich=3,
                    max_atomic_num=20,
                    save_path=save_dir,
                    oxidation_states_set=ox_states,
                )
                self.assertIsInstance(smact_df, pd.DataFrame)
                self.assertTrue(len(smact_df) == 1330)
                self.assertTrue(
                    smact_df["smact_allowed"].sum() == oxidation_states_sets_dict[ox_states]["smact_allowed"]
                )
                # Check if the data was saved to disk
                self.assertTrue(os.path.exists(save_dir))

                # Clean up
                shutil.rmtree("data")

    @pytest.mark.skipif(
        (
            sys.platform == "win32"
            or not (os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY"))
            or not MP_API_AVAILABLE
            or skip_mprester_tests
        ),
        reason="Test requires MP_API_KEY and fails on Windows due to filepath issues.",
    )
    def test_download_compounds_with_mp_api(self):
        save_mp_dir = "data/binary/mp_data"
        if MP_API_AVAILABLE:
            from smact.utils.crystal_space import download_compounds_with_mp_api

            download_compounds_with_mp_api.download_mp_data(
                mp_api_key=os.environ.get("MP_API_KEY"),
                num_elements=2,
                max_stoich=1,
                save_dir=save_mp_dir,
            )

        # Check if the data was downloaded
        self.assertTrue(os.path.exists(save_mp_dir))
        self.assertTrue(len(os.listdir(save_mp_dir)) > 0)

        # Clean up
        shutil.rmtree(save_mp_dir)


files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_ICSD_OX_STATES = os.path.join(files_dir, "oxidation_states_icsd24_consensus.txt")
TEST_ICSD_OX_STATES_W_ZERO = os.path.join(files_dir, "oxidation_states_icsd24_consensus_w_0.txt")


class OxidationStatesTest(unittest.TestCase):
    def setUp(self):
        self.ox_filter = ICSD24OxStatesFilter()
        with open(TEST_ICSD_OX_STATES) as f:
            self.test_ox_states = f.read()
        with open(TEST_ICSD_OX_STATES_W_ZERO) as f:
            self.test_ox_states_w_zero = f.read()

    def test_oxidation_states_filter(self):
        self.assertIsInstance(self.ox_filter.ox_states_df, pd.DataFrame)
        threshold = 10
        filtered_df = self.ox_filter.filter(consensus=threshold)

        self.assertIsInstance(filtered_df, pd.DataFrame)
        self.assertEqual(filtered_df.columns.tolist(), ["element", "oxidation_state"])

    def test_oxidation_states_write(self):
        self.maxDiff = None

        filename = "test_ox_states"
        comment = "Testing writing of ICSD 24 oxidation states list."
        self.ox_filter.write(
            filename,
            comment=comment,
            consensus=3,
            include_zero=False,
            commonality="low",
        )

        self.assertTrue(os.path.exists(f"{filename}.txt"))
        # Read the file and check its content
        with open(f"{filename}.txt") as f:
            content = f.read()
        # Check if the comment is included in the file
        self.assertIn(comment, content)
        # Check if the file content matches the expected content
        self.assertEqual(content, self.test_ox_states)
        os.remove(f"{filename}.txt")

        filename_w_zero = "test_ox_states_w_zero"
        comment = "Testing writing of ICSD 24 oxidation states list."
        self.ox_filter.write(
            filename_w_zero,
            comment=comment,
            consensus=3,
            include_zero=True,
            commonality="low",
        )

        self.assertTrue(os.path.exists(f"{filename_w_zero}.txt"))
        # Read the file and check its content
        with open(f"{filename_w_zero}.txt") as f:
            content = f.read()
        # Check if the comment is included in the file
        self.assertIn(comment, content)
        # Check if the file content matches the expected content
        self.assertEqual(content, self.test_ox_states_w_zero)
        os.remove(f"{filename_w_zero}.txt")

    def test_get_species_list(self):
        # Test with default parameters
        species_list = self.ox_filter.get_species_list()
        self.assertIsInstance(species_list, list)
        self.assertGreater(len(species_list), 0)  # Ensure the list is not empty

        # Test with include_zero=True
        species_list_with_zero = self.ox_filter.get_species_list(include_zero=True)
        self.assertIsInstance(species_list_with_zero, list)
        self.assertGreater(len(species_list_with_zero), 0)

        # Test with include_one_oxidation_state=True
        species_list_with_one = self.ox_filter.get_species_list(include_one_oxidation_state=True)
        self.assertIsInstance(species_list_with_one, list)
        self.assertGreater(len(species_list_with_one), 0)

        # Test with different commonality levels
        species_list_low = self.ox_filter.get_species_list(commonality="low")
        self.assertIsInstance(species_list_low, list)
        self.assertGreater(len(species_list_low), 0)

        species_list_medium = self.ox_filter.get_species_list(commonality="medium")
        self.assertIsInstance(species_list_medium, list)
        self.assertGreater(len(species_list_medium), 0)

        species_list_high = self.ox_filter.get_species_list(commonality="high")
        self.assertIsInstance(species_list_high, list)
        self.assertGreater(len(species_list_high), 0)

        # Test with a specific consensus threshold
        species_list_threshold = self.ox_filter.get_species_list(consensus=5)
        self.assertIsInstance(species_list_threshold, list)
        self.assertGreater(len(species_list_threshold), 0)

    def test_oxidation_states_filter_species_occurrences(self):
        species_occurrences_df = self.ox_filter.get_species_occurrences_df(consensus=1)
        self.assertIsInstance(species_occurrences_df, pd.DataFrame)
        self.assertEqual(
            species_occurrences_df.columns.tolist(),
            ["element", "species", "results_count", "species_proportion (%)"],
        )
        self.assertEqual(species_occurrences_df.shape, (490, 4))
        self.assertEqual(species_occurrences_df.iloc[0]["species"], "O2-")
        self.assertEqual(species_occurrences_df.iloc[0]["results_count"], 116910)
