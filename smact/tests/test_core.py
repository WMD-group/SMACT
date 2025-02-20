#!/usr/bin/env python
from __future__ import annotations

import os
import unittest

import pytest
from pymatgen.core import Composition, Structure
from pymatgen.core.periodic_table import Specie

import smact
import smact.distorter
import smact.lattice
import smact.lattice_parameters
import smact.oxidation_states
import smact.screening
from smact import Species
from smact.builder import wurtzite
from smact.properties import (
    band_gap_Harrison,
    compound_electroneg,
    valence_electron_count,
)

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_OX_STATES = os.path.join(files_dir, "test_oxidation_states.txt")
TEST_STRUCT = os.path.join(files_dir, "mp-540839_CsPbI3_oxi.json")


class TestSequenceFunctions(unittest.TestCase):
    # ---------------- TOP-LEVEL ----------------

    def test_Element_class_Pt(self):
        Pt = smact.Element(
            "Pt",
        )
        self.assertEqual(Pt.name, "Platinum")
        self.assertEqual(Pt.ionpot, 8.95883)
        self.assertEqual(Pt.number, 78)
        self.assertEqual(Pt.dipol, 44.00)

    def test_ordered_elements(self):
        self.assertEqual(smact.ordered_elements(65, 68), ["Tb", "Dy", "Ho", "Er"])
        self.assertEqual(smact.ordered_elements(52, 52), ["Te"])

    def test_element_dictionary(self):
        newlist = ["O", "Rb", "W"]
        dictionary = smact.element_dictionary(newlist, TEST_OX_STATES)
        self.assertEqual(dictionary["O"].crustal_abundance, 461000.0)
        self.assertEqual(dictionary["Rb"].oxidation_states_smact14, [-1, 1])
        self.assertEqual(dictionary["Rb"].oxidation_states, [1])
        self.assertEqual(dictionary["Rb"].oxidation_states_custom, [-1, 1])
        self.assertEqual(dictionary["W"].name, "Tungsten")
        self.assertTrue("Rn" in smact.element_dictionary())

    def test_are_eq(self):
        self.assertTrue(smact.are_eq([1.00, 2.00, 3.00], [1.001, 1.999, 3.00], tolerance=1e-2))
        self.assertFalse(smact.are_eq([1.00, 2.00, 3.00], [1.001, 1.999, 3.00]))

    def test_gcd_recursive(self):
        self.assertEqual(smact._gcd_recursive(4, 12, 10, 32), 2)
        self.assertEqual(smact._gcd_recursive(15, 4), 1)

    def test_isneutral(self):
        self.assertTrue(smact._isneutral((-2, 1), (2, 4)))
        self.assertFalse(smact._isneutral((4, -1), (1, 3)))

    def test_neutral_ratios(self):
        ox = [1, -2, 1]
        is_neutral, neutral_combos = smact.neutral_ratios(ox)
        self.assertTrue(is_neutral)
        self.assertEqual(len(neutral_combos), 9)
        self.assertTrue((3, 2, 1) in neutral_combos)

    # ---------------- Properties ----------------

    def test_compound_eneg_brass(self):
        self.assertAlmostEqual(
            compound_electroneg(elements=["Cu", "Zn"], stoichs=[0.5, 0.5], source="Pauling"),
            5.0638963259,
        )

    def test_harrison_gap_MgCl(self):
        self.assertAlmostEqual(
            band_gap_Harrison("Mg", "Cl", verbose=False, distance=2.67),
            3.545075110572662,
        )

    def test_valence_electron_count(self):
        # Test valid compounds
        self.assertAlmostEqual(valence_electron_count("Fe2O3"), 6.8, places=2)
        self.assertAlmostEqual(valence_electron_count("CuZn"), 11.5, places=2)

        # Test single element
        self.assertEqual(valence_electron_count("Fe"), 8)

        # Test empty string
        self.assertEqual(valence_electron_count(""), 0.0)

        # Test invalid elements and formats
        with pytest.raises(ValueError):
            valence_electron_count("Xx2O3")  # Xx is not a real element

        with pytest.raises(ValueError):
            valence_electron_count("LrO")

    # ---------------- BUILDER ----------------

    def test_builder_ZnS(self):
        ZnS, sys_ZnS = wurtzite(["Zn", "S"])
        self.assertEqual((ZnS.sites[0].position[2]), 0)
        self.assertEqual((ZnS.sites[0].position[0]), 2.0 / 3.0)

    # ---------------- SCREENING ----------------

    def test_pauling_test(self):
        """Test pauling_test functionality with various inputs and edge cases"""
        # Test basic functionality with valid inputs
        Sn, S = (smact.Element(lbl) for lbl in ("Sn", "S"))
        self.assertTrue(smact.screening.pauling_test([2, -2], [Sn.pauling_eneg, S.pauling_eneg], threshold=0.3))

        # Test with missing electronegativity
        self.assertFalse(smact.screening.pauling_test([2, -2], [None, S.pauling_eneg]))

        # Test with different oxidation states
        self.assertTrue(
            smact.screening.pauling_test(
                (+2, -2),
                (Sn.pauling_eneg, S.pauling_eneg),
            )
        )
        self.assertFalse(smact.screening.pauling_test((-2, +2), (Sn.pauling_eneg, S.pauling_eneg)))

        # Test with repeat anions/cations flags
        self.assertFalse(
            smact.screening.pauling_test(
                (-2, -2, +2),
                (S.pauling_eneg, S.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "S", "Sn"),
                repeat_anions=False,
            )
        )
        self.assertTrue(
            smact.screening.pauling_test(
                (-2, -2, +2),
                (S.pauling_eneg, S.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "S", "Sn"),
                repeat_cations=False,
            )
        )

    def test_pauling_test_old(self):
        Sn, S = (smact.Element(label) for label in ("Sn", "S"))
        self.assertTrue(
            smact.screening.pauling_test_old(
                (+2, -2),
                (Sn.pauling_eneg, S.pauling_eneg),
                symbols=("S", "S", "Sn"),
            )
        )
        self.assertFalse(
            smact.screening.pauling_test_old(
                (-2, +2),
                (Sn.pauling_eneg, S.pauling_eneg),
                symbols=("S", "S", "Sn"),
            )
        )
        self.assertFalse(
            smact.screening.pauling_test_old(
                (-2, -2, +2),
                (S.pauling_eneg, S.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "S", "Sn"),
                repeat_anions=False,
            )
        )
        self.assertTrue(
            smact.screening.pauling_test_old(
                (-2, -2, +2),
                (S.pauling_eneg, S.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "S", "Sn"),
                repeat_cations=False,
            )
        )
        self.assertFalse(
            smact.screening.pauling_test_old(
                (-2, +2, +2),
                (S.pauling_eneg, Sn.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "Sn", "Sn"),
                repeat_cations=False,
            )
        )
        self.assertTrue(
            smact.screening.pauling_test_old(
                (-2, +2, +2),
                (S.pauling_eneg, Sn.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "Sn", "Sn"),
                repeat_anions=False,
            )
        )

    def test_eneg_states_test(self):
        Na, Fe, Cl = (smact.Element(label) for label in ("Na", "Fe", "Cl"))
        self.assertTrue(
            smact.screening.eneg_states_test([1, 3, -1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg])
        )
        self.assertFalse(
            smact.screening.eneg_states_test([-1, 3, 1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg])
        )

    def test_eneg_states_test_alternate(self):
        Na, Fe, Cl = (smact.Element(label) for label in ("Na", "Fe", "Cl"))
        self.assertTrue(
            smact.screening.eneg_states_test_alternate([1, 3, -1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg])
        )
        self.assertFalse(
            smact.screening.eneg_states_test_alternate([-1, 3, 1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg])
        )

    def test_eneg_states_test_threshold(self):
        self.assertFalse(smact.screening.eneg_states_test_threshold([1, -1], [1.83, 1.82], threshold=0))
        self.assertTrue(smact.screening.eneg_states_test_threshold([1, -1], [1.83, 1.82], threshold=0.1))

    def test_ml_rep_generator(self):
        Pb, O = (smact.Element(label) for label in ("Pb", "O"))
        PbO2_ml = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.6666666666666666,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.3333333333333333,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        self.assertEqual(smact.screening.ml_rep_generator(["Pb", "O"], [1, 2]), PbO2_ml)
        self.assertEqual(smact.screening.ml_rep_generator([Pb, O], [1, 2]), PbO2_ml)

    def test_smact_filter(self):
        oxidation_states_sets = ["smact14", "icsd24"]
        oxidation_states_sets_results = {
            "smact14": {
                "thresh_2": [
                    (("Na", "Fe", "Cl"), (1, -1, -1), (2, 1, 1)),
                    (("Na", "Fe", "Cl"), (1, 1, -1), (1, 1, 2)),
                ]
            },
            "icsd24": {"thresh_2": [(("Na", "Fe", "Cl"), (1, 1, -1), (1, 1, 2))]},
        }

        Na, Fe, Cl = (smact.Element(label) for label in ("Na", "Fe", "Cl"))

        for ox_state_set in oxidation_states_sets:
            with self.subTest(ox_state_set=ox_state_set):
                output = smact.screening.smact_filter([Na, Fe, Cl], threshold=2, oxidation_states_set=ox_state_set)
                self.assertEqual(
                    [(r[0], r[1], r[2]) for r in output],
                    oxidation_states_sets_results[ox_state_set]["thresh_2"],
                )

        # Test that reading the oxidation states from a file produces the same results
        self.assertEqual(
            smact.screening.smact_filter([Na, Fe, Cl], threshold=2, oxidation_states_set="smact14"),
            smact.screening.smact_filter([Na, Fe, Cl], threshold=2, oxidation_states_set=TEST_OX_STATES),
        )

        self.assertEqual(
            set(
                smact.screening.smact_filter(
                    [Na, Fe, Cl], threshold=2, species_unique=False, oxidation_states_set="smact14"
                )
            ),
            {
                (("Na", "Fe", "Cl"), (2, 1, 1)),
                (("Na", "Fe", "Cl"), (1, 1, 2)),
            },
        )

        self.assertEqual(
            len(smact.screening.smact_filter([Na, Fe, Cl], threshold=8, oxidation_states_set="smact14")), 77
        )

        result = smact.screening.smact_filter([Na, Fe, Cl], stoichs=[[1], [1], [4]], oxidation_states_set="smact14")
        self.assertEqual(
            [(r[0], r[1], r[2]) for r in result],
            [
                (("Na", "Fe", "Cl"), (1, 3, -1), (1, 1, 4)),
            ],
        )
        stoichs = [list(range(1, 5)), list(range(1, 5)), list(range(1, 10))]
        self.assertEqual(
            len(smact.screening.smact_filter([Na, Fe, Cl], stoichs=stoichs, oxidation_states_set="smact14")),
            45,
        )

    # --------- Tests for Screening Logic ---------

    def test_get_valid_combinations(self):
        """Test get_valid_combinations returns valid (element_symbols, ratio) combinations"""
        valid_combos = smact.screening.get_valid_combinations("NaCl")
        self.assertIsInstance(valid_combos, list)
        self.assertGreater(len(valid_combos), 0)
        for combo in valid_combos:
            self.assertEqual(len(combo), 2)
            self.assertIsInstance(combo[0], tuple)
            self.assertIsInstance(combo[1], tuple)

    def test_smact_validity_return_all(self):
        """Test smact_validity with return_all=True returns valid combinations"""
        result = smact.screening.smact_validity("NaCl", return_all=True)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)
        for combo in result:
            self.assertEqual(len(combo), 2)
            self.assertIsInstance(combo[0], tuple)
            self.assertIsInstance(combo[1], tuple)

    def test_smact_validity_short_circuit(self):
        """Test smact_validity with return_all=False returns boolean"""
        result = smact.screening.smact_validity("NaCl", return_all=False)
        self.assertIsInstance(result, bool)
        self.assertTrue(result)

    def test_smact_validity_invalid(self):
        """Test invalid composition returns False or empty list"""
        result_bool = smact.screening.smact_validity("Al3Li", include_alloys=False, return_all=False)
        self.assertFalse(result_bool)
        result_all = smact.screening.smact_validity("Al3Li", include_alloys=False, return_all=True)
        self.assertEqual(result_all, [])

    def test_generate_valid_combos(self):
        """Test _generate_valid_combos yields valid tuples"""
        from pymatgen.core import Composition

        combos = list(smact.screening._generate_valid_combos(Composition("NaCl")))
        self.assertGreater(len(combos), 0)
        for combo in combos:
            self.assertEqual(len(combo), 3)
            self.assertIsInstance(combo[0], tuple)
            self.assertIsInstance(combo[1], tuple)
            self.assertIsInstance(combo[2], tuple)

        # Test additional branches
        single_elem_combos = list(smact.screening._generate_valid_combos(Composition("Fe")))
        self.assertEqual(len(single_elem_combos), 1)
        self.assertEqual(single_elem_combos[0], (("Fe",), (0,), (1,)))

        alloy_combos = list(smact.screening._generate_valid_combos(Composition("FeAl"), include_alloys=True))
        self.assertGreater(len(alloy_combos), 0)

        metallic_combos = list(
            smact.screening._generate_valid_combos(Composition("Fe"), check_metallicity=True, metallicity_threshold=0.5)
        )
        self.assertGreater(len(metallic_combos), 0)

        no_pauling_combos = list(smact.screening._generate_valid_combos(Composition("NaCl"), use_pauling_test=False))
        self.assertGreater(len(no_pauling_combos), 0)

    def test_pauling_test_type_error_handling(self):
        """Test handling of TypeError in pauling test"""
        import smact.screening

        def mock_pauling_test(*args, **kwargs):
            raise TypeError("Mock TypeError")

        original_pauling_test = smact.screening.pauling_test
        try:
            smact.screening.pauling_test = mock_pauling_test
            result = smact.screening.smact_validity("NaCl", use_pauling_test=True, return_all=False)
            self.assertTrue(result)
        finally:
            smact.screening.pauling_test = original_pauling_test

    def test_smact_validity_type_error_electronegativity(self):
        """Test pauling_test handles missing electronegativity"""
        He = smact.Element("He")
        original_eneg = He.pauling_eneg
        He.pauling_eneg = None
        try:
            result = smact.screening.smact_validity("HeF2")
            self.assertIsInstance(result, bool)
        finally:
            He.pauling_eneg = original_eneg

    def test_type_checking(self):
        """Test type checking and annotations"""
        from collections.abc import Generator
        from typing import get_type_hints

        # Test return type hints
        hints = get_type_hints(smact.screening._generate_valid_combos)
        self.assertIn("return", hints)
        self.assertTrue(issubclass(hints["return"], Generator))

        # Test generator instance
        combos = smact.screening._generate_valid_combos(Composition("Fe"))
        self.assertIsInstance(combos, Generator)
        combo_list = list(combos)
        self.assertGreater(len(combo_list), 0)

    # ---------------- Lattice ----------------
    def test_Lattice_class(self):
        site_A = smact.lattice.Site([0, 0, 0], -1)
        site_B = smact.lattice.Site([0.5, 0.5, 0.5], [+2, +3])
        test_lattice = smact.lattice.Lattice([site_A, site_B], space_group=221)
        self.assertEqual(test_lattice.sites[0].oxidation_states, -1)
        self.assertEqual(test_lattice.sites[1].position, [0.5, 0.5, 0.5])

    # ---------- Lattice parameters -----------
    def test_lattice_parameters(self):
        perovskite = smact.lattice_parameters.cubic_perovskite([1.81, 1.33, 1.82])
        wurtz = smact.lattice_parameters.wurtzite([1.81, 1.33])
        self.assertAlmostEqual(perovskite[0], 6.3)
        self.assertAlmostEqual(perovskite[1], 6.3)
        self.assertAlmostEqual(perovskite[3], 90)
        self.assertAlmostEqual(wurtz[0], 5.13076)
        self.assertAlmostEqual(wurtz[2], 8.3838)

    # ---------- smact.oxidation_states module -----------
    def test_oxidation_states(self):
        ox = smact.oxidation_states.Oxidation_state_probability_finder()
        self.assertAlmostEqual(
            ox.compound_probability([Specie("Fe", +3), Specie("O", -2)]),
            0.74280230326,
        )
        self.assertAlmostEqual(
            ox.pair_probability(Species("Fe", +3), Species("O", -2)),
            0.74280230326,
        )
        self.assertEqual(len(ox.get_included_species()), 173)

    def test_compound_probability_structure(self):
        structure = Structure.from_file(TEST_STRUCT)
        ox = smact.oxidation_states.Oxidation_state_probability_finder()
        self.assertEqual(ox.compound_probability(structure), 1.0)
