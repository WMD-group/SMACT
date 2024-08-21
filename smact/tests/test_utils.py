import os
import unittest

import pandas as pd
import pytest

from smact.utils.oxidation import ICSD24OxStatesFilter

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_ICSD_OX_STATES = os.path.join(
    files_dir, "test_icsd_oxidation_states_filter_1000.txt"
)
TEST_ICSD_OX_STATES_W_ZERO = os.path.join(
    files_dir, "test_icsd_oxidation_states_filter_1000_w_0_ox_state.txt"
)


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
        filtered_df = self.ox_filter.filter(threshold)

        self.assertIsInstance(filtered_df, pd.DataFrame)
        self.assertEqual(
            filtered_df.columns.tolist(), ["element", "oxidation_state"]
        )
        # self.assertEqual(filtered_df.loc[""])

    def test_oxidation_states_write(self):
        threshold = 1000
        filename = "test_ox_states"
        filename_w_zero = "test_ox_states_w_zero"
        comment = "Testing writing of ICSD 24 oxidation states list."
        self.ox_filter.write(filename, threshold, comment=comment)
        self.ox_filter.write(
            filename_w_zero, threshold, include_zero=True, comment=comment
        )
        self.assertTrue(os.path.exists(f"{filename}.txt"))
        assert [line for line in open(f"{filename}.txt")] == [
            line for line in open(TEST_ICSD_OX_STATES)
        ]
        self.assertTrue(os.path.exists(f"{filename}_w_zero.txt"))
        assert [line for line in open(f"{filename}_w_zero.txt")] == [
            line for line in open(TEST_ICSD_OX_STATES_W_ZERO)
        ]
        # Clean up
        os.remove(f"{filename}.txt")
        os.remove(f"{filename}_w_zero.txt")

    def test_oxidation_states_filter_species_list(self):
        for threshold, length in [(0, 460), (5, 337), (50, 214)]:
            species_list = self.ox_filter.get_species_list(threshold)
            self.assertIsInstance(species_list, list)
            self.assertEqual(len(species_list), length)
