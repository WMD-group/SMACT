import os
import unittest

import smact
from smact.dopant_prediction import doper
from smact.structure_prediction import mutation, utilities

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_LAMBDA_JSON = os.path.join(files_dir, "test_lambda_tab.json")


class DopantPredictionTest(unittest.TestCase):
    def test_dopant_prediction(self):
        num_dopants = 10
        test_specie = ("Cu+", "Ga3+", "S2-")
        test = doper.Doper(test_specie)

        cation_max_charge = max(
            test_specie, key=lambda x: utilities.parse_spec(x)[1]
        )
        anion_min_charge = min(
            test_specie, key=lambda x: utilities.parse_spec(x)[1]
        )
        _, cat_charge = utilities.parse_spec(cation_max_charge)
        _, an_charge = utilities.parse_spec(anion_min_charge)

        # Assert: Length of the list and return type (dictionary: list)
        dopants = test.get_dopants()
        result = test.results
        self.assertIs(type(result), dict)
        for d in result.values():
            self.assertIn("sorted", d)
            for v in d.values():
                self.assertIs(type(v), list)

        # Assert: (cation) higher charges for n-type and lower charges for p-type
        n_sub_list_cat = result.get("n-type cation substitutions").get("sorted")
        p_sub_list_cat = result.get("p-type cation substitutions").get("sorted")
        n_sub_list_an = result.get("n-type anion substitutions").get("sorted")
        p_sub_list_an = result.get("p-type anion substitutions").get("sorted")

        for n_atom, p_atom in zip(n_sub_list_cat, p_sub_list_cat):
            self.assertGreater(utilities.parse_spec(n_atom[0])[1], cat_charge)
            self.assertLess(utilities.parse_spec(p_atom[0])[1], cat_charge)

        for n_atom, p_atom in zip(n_sub_list_an, p_sub_list_an):
            self.assertGreater(utilities.parse_spec(n_atom[0])[1], an_charge)
            self.assertLess(utilities.parse_spec(p_atom[0])[1], an_charge)

    def test_dopant_prediction_skipspecies(self):
        test_specie = ("Cu+", "Ga3+", "S2-")
        with self.assertRaises(ValueError):
            doper.Doper(
                test_specie, filepath=TEST_LAMBDA_JSON, embedding="skipspecies"
            )

        with self.assertRaises(ValueError):
            doper.Doper(test_specie, embedding="skip", use_probability=False)

        test = doper.Doper(
            test_specie, embedding="skipspecies", use_probability=False
        )
        dopants = test.get_dopants()
        result = test.results

        n_sub_list_cat = result.get("n-type cation substitutions").get("sorted")
        p_sub_list_cat = result.get("p-type cation substitutions").get("sorted")
        n_sub_list_an = result.get("n-type anion substitutions").get("sorted")
        p_sub_list_an = result.get("p-type anion substitutions").get("sorted")

        # Create a list of the results
        results = [n_sub_list_cat, p_sub_list_cat, n_sub_list_an, p_sub_list_an]

        # Assert that the similarity values are between 0 and 1
        for r in results:
            for dopant_result_list in r:
                self.assertTrue(0 <= dopant_result_list[2] <= 1)

    def test_format_number(self):
        test_specie = ("Cu+", "Ga3+", "S2-")
        test = doper.Doper(test_specie)

        self.assertEqual(test._format_number(2), "2+")
        self.assertEqual(test._format_number(-2), "2-")
