from __future__ import annotations

import os
import shutil
import sys
import unittest
from importlib.util import find_spec
from unittest.mock import patch

import pandas as pd
import pytest
import requests
from pymatgen.core import SETTINGS, Composition

from smact import Element
from smact.data_loader import (
    lookup_element_oxidation_states_custom,
    lookup_element_shannon_radius_data_extendedML,
    lookup_element_sse2015_data,
)
from smact.screening import SmactFilterOutputs, smact_filter
from smact.structure_prediction.database import _validate_table_name
from smact.utils.composition import comp_maker, composition_dict_maker, formula_maker, parse_formula
from smact.utils.crystal_space import download_compounds_with_mp_api, generate_composition_with_smact
from smact.utils.oxidation import ICSD24OxStatesFilter
from smact.utils.species import _parse_spec_old, parse_spec

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
        comp4 = comp_maker(self.smact_filter_output[1])  # type: ignore[arg-type]
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
        form4 = formula_maker(self.smact_filter_output[1])  # type: ignore[arg-type]
        self.assertEqual(form1, "FeO")
        self.assertEqual(form2, "FeO")
        self.assertEqual(form1, form2)
        self.assertEqual(form3, "Fe3O4")
        self.assertEqual(form4, "Li10Ge(PS6)2")

    def test_composition_dict_maker(self):
        """Test the composition_dict_maker function"""
        # 3-element tuple (symbols, ox_states, stoichs) - returns species keys
        dict1 = composition_dict_maker(self.mock_filter_output[0])
        self.assertIsInstance(dict1, dict)
        self.assertIn("Fe2+", dict1)
        self.assertIn("O2-", dict1)
        self.assertEqual(dict1["Fe2+"], 1.0)
        self.assertEqual(dict1["O2-"], 1.0)

        # 2-element tuple (symbols, stoichs) - returns element keys
        dict2 = composition_dict_maker(self.mock_filter_output[1])
        self.assertIsInstance(dict2, dict)
        self.assertIn("Fe", dict2)
        self.assertIn("O", dict2)

    def test_smact_filter_return_output_default(self):
        """Test smact_filter with default return_output returns tuples"""
        els = [Element("Li"), Element("O")]
        comps = smact_filter(els, threshold=2)
        self.assertIsInstance(comps, list)
        self.assertGreater(len(comps), 0)
        # Default should return tuples
        self.assertIsInstance(comps[0], tuple)

    def test_smact_filter_return_output_formula(self):
        """Test smact_filter with return_output='formula' returns strings"""
        els = [Element("Li"), Element("O")]
        comps = smact_filter(els, threshold=2, return_output=SmactFilterOutputs.formula)
        self.assertIsInstance(comps, list)
        self.assertGreater(len(comps), 0)
        for comp in comps:
            self.assertIsInstance(comp, str)
        self.assertIn("Li2O", comps)

    def test_smact_filter_return_output_dict(self):
        """Test smact_filter with return_output='dict' returns dicts"""
        els = [Element("Li"), Element("O")]
        comps = smact_filter(els, threshold=2, return_output=SmactFilterOutputs.composition_dict)
        self.assertIsInstance(comps, list)
        self.assertGreater(len(comps), 0)
        for comp in comps:
            self.assertIsInstance(comp, dict)

    def test_smact_filter_return_output_species_not_unique(self):
        """Test smact_filter return_output with species_unique=False"""
        els = [Element("Li"), Element("O")]
        # Formula output
        formulas = smact_filter(
            els,
            threshold=2,
            species_unique=False,
            return_output=SmactFilterOutputs.formula,
        )
        self.assertIsInstance(formulas, list)
        self.assertGreater(len(formulas), 0)
        for f in formulas:
            self.assertIsInstance(f, str)

        # Dict output
        dicts = smact_filter(
            els,
            threshold=2,
            species_unique=False,
            return_output=SmactFilterOutputs.composition_dict,
        )
        self.assertIsInstance(dicts, list)
        self.assertGreater(len(dicts), 0)
        for d in dicts:
            self.assertIsInstance(d, dict)

    def test_smact_filter_invalid_return_output(self):
        """Test smact_filter raises ValueError for invalid return_output"""
        els = [Element("Li"), Element("O")]
        with pytest.raises(ValueError):
            smact_filter(els, threshold=2, return_output="invalid")  # type: ignore[arg-type]


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

    def test_generate_composition_with_smact_custom(self):
        save_dir = "data/binary/df_binary_label.pkl"
        oxidation_states_sets = ["smact14", "icsd24"]
        oxidation_states_sets_dict = {
            "smact14": {"smact_allowed": 388},
            "icsd24": {"smact_allowed": 342},
        }
        for ox_states in oxidation_states_sets:
            with self.subTest(ox_states=ox_states):
                smact_df = generate_composition_with_smact.generate_composition_with_smact_custom(
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

        # Test commonality="main" (should return species with max proportion for each element)
        species_list_main = self.ox_filter.get_species_list(commonality="main")
        self.assertIsInstance(species_list_main, list)
        self.assertGreater(len(species_list_main), 0)

        # Get the original dataframe for comparison
        df_main = self.ox_filter.get_species_occurrences_df()
        max_proportions = df_main.groupby("element")["species_proportion (%)"].max()  # type: ignore[union-attr]

        # Verify that only species with maximum proportion for each element are included
        for species in species_list_main:
            element = species.split("+")[0].split("-")[0].rstrip("0123456789")  # Extract element from species string
            species_proportion = df_main[df_main["species"] == species]["species_proportion (%)"].iloc[0]  # type: ignore[union-attr]
            self.assertEqual(species_proportion, max_proportions[element])

        # Test specific commonality threshold
        threshold = 50.0  # Testing with 50% threshold
        species_list_threshold = self.ox_filter.get_species_list(commonality=threshold)
        self.assertIsInstance(species_list_threshold, list)

        # Verify that all species meet the threshold requirement
        df_threshold = self.ox_filter.get_species_occurrences_df()
        for species in species_list_threshold:
            proportion = df_threshold[df_threshold["species"] == species]["species_proportion (%)"].iloc[0]  # type: ignore[union-attr]
            self.assertGreaterEqual(proportion, threshold)

        # Test that "main" returns different results than threshold-based filtering
        species_list_main = self.ox_filter.get_species_list(commonality="main")
        species_list_high_threshold = self.ox_filter.get_species_list(commonality=90.0)
        self.assertNotEqual(set(species_list_main), set(species_list_high_threshold))

    def test_oxidation_states_filter_species_occurrences(self):
        species_occurrences_df = self.ox_filter.get_species_occurrences_df(consensus=1)
        self.assertIsInstance(species_occurrences_df, pd.DataFrame)
        self.assertEqual(
            species_occurrences_df.columns.tolist(),  # type: ignore[union-attr]
            ["element", "species", "results_count", "species_proportion (%)"],
        )
        self.assertEqual(species_occurrences_df.shape, (490, 4))
        self.assertEqual(species_occurrences_df.iloc[0]["species"], "O2-")  # type: ignore[union-attr]
        self.assertEqual(species_occurrences_df.iloc[0]["results_count"], 116910)  # type: ignore[union-attr]


class TestCompositionEdgeCases(unittest.TestCase):
    """Edge case coverage for smact/utils/composition.py."""

    def test_parse_formula_invalid_chars(self):
        """_get_sym_dict raises ValueError for invalid formula (line 49-50)."""
        with pytest.raises(ValueError, match="invalid formula"):
            parse_formula("123invalid")

    def test_parse_formula_empty_string(self):
        """parse_formula with empty string returns empty dict."""
        result = parse_formula("")
        self.assertEqual(result, {})


class TestOxidationEdgeCases(unittest.TestCase):
    """Edge case coverage for ICSD24OxStatesFilter."""

    def setUp(self):
        self.ox_filter = ICSD24OxStatesFilter()

    def test_commonality_type_error(self):
        """filter raises TypeError for invalid commonality type (line 53)."""
        with pytest.raises(TypeError, match="commonality must be"):
            self.ox_filter.filter(commonality=[1, 2, 3])  # type: ignore[arg-type]

    def test_get_species_occurrences_unsorted(self):
        """get_species_occurrences_df with sort_by_occurrences=False (line 157)."""
        df = self.ox_filter.get_species_occurrences_df(sort_by_occurrences=False)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)  # type: ignore[arg-type]

    def test_filter_numeric_commonality(self):
        """filter with a numeric commonality threshold (line 51)."""
        df = self.ox_filter.filter(commonality=25.0)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)

    def test_get_species_occurrences_include_zero(self):
        """get_species_occurrences_df with include_zero=True (line 138)."""
        df = self.ox_filter.get_species_occurrences_df(include_zero=True)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)  # type: ignore[arg-type]

    def test_get_species_list_value_error_path(self):
        """get_species_list skips species whose ox_state can't be parsed (lines 113-114)."""
        # Inject a corrupted row where oxidation_state is a non-integer string.
        # The code splits on spaces, so use a single token that fails int().
        corrupted_df = pd.DataFrame(
            {
                "element": ["Na", "Na"],
                "oxidation_state": ["bad", "1"],
            }
        )
        with patch.object(self.ox_filter, "filter", return_value=corrupted_df):
            species_list = self.ox_filter.get_species_list()

        # int("bad") → ValueError → skipped via continue
        # Only "Na1+" from the second row should appear
        self.assertEqual(len(species_list), 1)
        self.assertIn("Na+", species_list)


class TestDatabaseEdgeCases(unittest.TestCase):
    """Edge case coverage for database.py _validate_table_name."""

    def test_invalid_table_name(self):
        """_validate_table_name raises ValueError for SQL-unsafe names (line 46)."""
        with pytest.raises(ValueError, match="Invalid table name"):
            _validate_table_name("DROP TABLE; --")


class TestSpeciesParsing(unittest.TestCase):
    """Branch coverage for smact/utils/species.py _parse_spec_old."""

    def setUp(self):
        self._parse_spec_old = _parse_spec_old

    def test_parse_spec_old_no_element_symbol_raises(self):
        """Line 38: string starting with digit has no element symbol → ValueError."""
        with pytest.raises(ValueError, match="Invalid species string"):
            self._parse_spec_old("123")

    def test_parse_spec_old_negative_ox_state(self):
        """Line 45: '-' in species with a digit → ox_state *= -1."""
        # "Fe(-2)": main regex fails (parens), _parse_spec_old finds digit '2' and '-'
        ele, charge = self._parse_spec_old("Fe(-2)")
        self.assertEqual(ele, "Fe")
        self.assertEqual(charge, -2)

    def test_parse_spec_old_zero_with_digit_zero(self):
        """Line 51: ox_state==0 and '0' in species → stays 0."""
        # "Fe0": no sign → main regex fails → _parse_spec_old
        ele, charge = parse_spec("Fe0")
        self.assertEqual(ele, "Fe")
        self.assertEqual(charge, 0)

    def test_parse_spec_old_bare_plus(self):
        """Line 53: '+' in species with no digit → ox_state=1."""
        # "Fe(+)": parens block main regex match; no digit; '+' present
        ele, charge = self._parse_spec_old("Fe(+)")
        self.assertEqual(ele, "Fe")
        self.assertEqual(charge, 1)

    def test_parse_spec_old_bare_minus(self):
        """Line 55: '-' in species, no digit, no '+', no '0' → ox_state=-1."""
        # "Fe(-)": parens block main regex; no digit; '-' present
        ele, charge = self._parse_spec_old("Fe(-)")
        self.assertEqual(ele, "Fe")
        self.assertEqual(charge, -1)


class TestDataLoaderWarnings(unittest.TestCase):
    """Branch coverage for smact/data_loader.py warning and missing-symbol paths."""

    def test_warn_on_missing_logs_debug(self):
        """_warn logs a debug message for missing element lookups."""
        test_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files", "test_oxidation_states.txt")
        with self.assertLogs("smact.data_loader", level="DEBUG") as cm:
            result = lookup_element_oxidation_states_custom("Xx", test_file)
        self.assertIsNone(result)
        self.assertTrue(any("not found" in msg for msg in cm.output))

    def test_shannon_radii_extendedML_found_and_missing(self):
        """Lines 381-398 (loader body), 447-451 (missing symbol path)."""
        # Valid element: exercises _load_shannon_radii_extendedML (lines 381-398)
        data = lookup_element_shannon_radius_data_extendedML("Fe")
        self.assertIsNotNone(data)
        self.assertIsInstance(data, list)

        # Unknown element: exercises missing-symbol warning path (lines 447-451)
        result = lookup_element_shannon_radius_data_extendedML("Xx")
        self.assertIsNone(result)

    def test_sse2015_missing_symbol(self):
        """Lines 560-561: lookup_element_sse2015_data with unknown symbol returns None."""
        result = lookup_element_sse2015_data("Xx")
        self.assertIsNone(result)
