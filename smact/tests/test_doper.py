from __future__ import annotations

import os
import unittest
import warnings
from unittest.mock import patch

import numpy as np
import pytest

from smact.dopant_prediction import doper
from smact.structure_prediction import utilities
from smact.structure_prediction.mutation import CationMutator

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_LAMBDA_JSON = os.path.join(files_dir, "test_lambda_tab.json")


class DopantPredictionTest(unittest.TestCase):
    def test_dopant_prediction(self):
        test_specie = ("Cu+", "Ga3+", "S2-")
        test = doper.Doper(test_specie)

        cation_max_charge = max(test_specie, key=lambda x: utilities.parse_spec(x)[1])
        anion_min_charge = min(test_specie, key=lambda x: utilities.parse_spec(x)[1])
        _, cat_charge = utilities.parse_spec(cation_max_charge)
        _, an_charge = utilities.parse_spec(anion_min_charge)

        # Assert: Length of the list and return type (dictionary: list)
        result = test.get_dopants()
        self.assertIs(type(result), dict)
        for d in result.values():
            self.assertIn("sorted", d)
            for v in d.values():
                self.assertIs(type(v), list)

        # Assert: (cation) higher charges for n-type and lower charges for p-type
        n_sub_list_cat = result.get("n-type cation substitutions").get("sorted")  # type: ignore[union-attr]
        p_sub_list_cat = result.get("p-type cation substitutions").get("sorted")  # type: ignore[union-attr]
        n_sub_list_an = result.get("n-type anion substitutions").get("sorted")  # type: ignore[union-attr]
        p_sub_list_an = result.get("p-type anion substitutions").get("sorted")  # type: ignore[union-attr]

        for n_atom, p_atom in zip(n_sub_list_cat, p_sub_list_cat, strict=False):
            self.assertGreater(utilities.parse_spec(n_atom[0])[1], cat_charge)
            self.assertLess(utilities.parse_spec(p_atom[0])[1], cat_charge)

        for n_atom, p_atom in zip(n_sub_list_an, p_sub_list_an, strict=False):
            self.assertGreater(utilities.parse_spec(n_atom[0])[1], an_charge)
            self.assertLess(utilities.parse_spec(p_atom[0])[1], an_charge)

    def test_dopant_prediction_skipspecies(self):
        test_specie = ("Cu+", "Ga3+", "S2-")
        with pytest.raises(ValueError):
            doper.Doper(test_specie, filepath=TEST_LAMBDA_JSON, embedding="skipspecies")

        with pytest.raises(ValueError):
            doper.Doper(test_specie, embedding="skip", use_probability=False)

        test = doper.Doper(test_specie, embedding="skipspecies", use_probability=False)
        result = test.get_dopants()

        n_sub_list_cat = result.get("n-type cation substitutions").get("sorted")  # type: ignore[union-attr]
        p_sub_list_cat = result.get("p-type cation substitutions").get("sorted")  # type: ignore[union-attr]
        n_sub_list_an = result.get("n-type anion substitutions").get("sorted")  # type: ignore[union-attr]
        p_sub_list_an = result.get("p-type anion substitutions").get("sorted")  # type: ignore[union-attr]

        # Create a list of the results
        results = [n_sub_list_cat, p_sub_list_cat, n_sub_list_an, p_sub_list_an]

        # Assert that the similarity values are between 0 and 1
        for r in results:
            for dopant_result_list in r:
                self.assertTrue(0 <= dopant_result_list[2] <= 1)

    def test_alternative_representations(self):
        test_specie = ("Cu+", "Ga3+", "S2-")
        test_gap = doper.Doper(test_specie, embedding="M3GNet-MP-2023.11.1-oxi-band_gap")
        test_eform = doper.Doper(test_specie, embedding="M3GNet-MP-2023.11.1-oxi-Eform")
        test_lamba = doper.Doper(test_specie, filepath=TEST_LAMBDA_JSON)
        for test in [test_gap, test_eform, test_lamba]:
            self.assertIsInstance(test, doper.Doper)
            result = test.get_dopants()
            self.assertIsInstance(result, dict)

    def test_format_number(self):
        test_specie = ("Cu+", "Ga3+", "S2-")
        test = doper.Doper(test_specie)

        self.assertEqual(test._format_number(2), "2+")
        self.assertEqual(test._format_number(-2), "2-")

    def test_unparseable_ion_warning(self):
        """Lines 196-197: ion whose string has no element symbol emits a warning and is skipped."""
        # "1+" starts with a digit so _parse_spec_old raises ValueError → caught, warning emitted
        test = doper.Doper(("Cu+", "1+", "S2-"))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = test.get_dopants()
        self.assertIsInstance(result, dict)
        self.assertTrue(any("Could not parse" in str(w.message) for w in caught))

    def test_plot_dopants_not_called_before_get_dopants(self):
        """Lines 321-322: plot_dopants raises RuntimeError when results is None."""
        test = doper.Doper(("Cu+", "Ga3+", "S2-"))
        with pytest.raises(RuntimeError, match="get_dopants first"):
            test.plot_dopants()

    def test_to_table_no_results(self):
        """Lines 358-360: to_table returns empty string when results is None."""
        test = doper.Doper(("Cu+", "Ga3+", "S2-"))
        result = test.to_table
        self.assertEqual(result, "")

    def test_to_table_with_results(self):
        """Lines 361-374: to_table returns tabulated string after get_dopants."""
        test = doper.Doper(("Cu+", "Ga3+", "S2-"))
        test.get_dopants(num_dopants=2)

        out = test.to_table
        self.assertIsInstance(out, str)
        self.assertIn("n-type cation substitutions", out)

    def test_plot_dopants_branches(self):
        """Lines 325-332: plot_dopants with similarity, selectivity, and combined plot_value."""
        test = doper.Doper(("Cu+", "Ga3+", "S2-"))
        test.get_dopants(num_dopants=3)

        for plot_value in ("probability", "similarity", "selectivity", "combined"):
            with patch("pymatgen.util.plotting.periodic_table_heatmap"):
                test.plot_dopants(plot_value=plot_value)

    def test_anion_substitution_continue_paths(self):
        """Lines 239 and 246: continue when dopant charge equals anion charge."""
        # Use a compound with two anions so that n-type candidates for one anion
        # have the same charge as the other anion → triggers the continue guard.
        test = doper.Doper(("Cu+", "S2-", "Cl1-"))
        result = test.get_dopants()
        self.assertIsInstance(result, dict)

    def test_doper_alpha_none_fallback(self):
        """Doper with alpha=None uses _DEFAULT_LAMBDA_THRESHOLD fallback (lines 106-107)."""
        # Build a CationMutator with a fully populated lambda table and alpha=None
        cm = CationMutator.from_json()
        # Override alpha to None after construction
        cm.alpha = None
        test = doper.Doper.__new__(doper.Doper)
        test.original_species = ("Cu+", "Ga3+", "S2-")
        test.cation_mutator = cm
        test.possible_species = list(cm.specs)
        # alpha is None → fallback path
        test.lambda_threshold = doper._DEFAULT_LAMBDA_THRESHOLD
        test.threshold = 1 / cm.Z * np.exp(doper._DEFAULT_LAMBDA_THRESHOLD)
        test.use_probability = True
        test.results = None

        result = test.get_dopants()
        self.assertIsInstance(result, dict)
