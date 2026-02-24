#!/usr/bin/env python
from __future__ import annotations

import os
import unittest
from os.path import dirname, exists, join, realpath

import pytest
from pymatgen.core import Structure
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
        ZnS, _ = wurtzite(["Zn", "S"])
        self.assertEqual((ZnS.sites[0].position[2]), 0)
        self.assertEqual((ZnS.sites[0].position[0]), 2.0 / 3.0)

    # ---------------- SCREENING ----------------

    def test_pauling_test(self):
        Sn, S = (smact.Element(label) for label in ("Sn", "S"))
        self.assertTrue(
            smact.screening.pauling_test(
                (+2, -2),
                (Sn.pauling_eneg, S.pauling_eneg),
            )
        )
        self.assertFalse(smact.screening.pauling_test((-2, +2), (Sn.pauling_eneg, S.pauling_eneg)))
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
        self.assertFalse(
            smact.screening.pauling_test(
                (-2, +2, +2),
                (S.pauling_eneg, Sn.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "Sn", "Sn"),
                repeat_cations=False,
            )
        )
        self.assertTrue(
            smact.screening.pauling_test(
                (-2, +2, +2),
                (S.pauling_eneg, Sn.pauling_eneg, Sn.pauling_eneg),
                symbols=("S", "Sn", "Sn"),
                repeat_anions=False,
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

    # --------- New tests for revised screening logic ---------

    def test_smact_validity(self):
        """
        Test that smact_validity returns a boolean indicating whether a composition is valid.
        """
        # Test with a valid composition
        result = smact.screening.smact_validity("NaCl")
        self.assertIsInstance(result, bool)
        self.assertTrue(result)

        # Test with an invalid composition
        invalid_result = smact.screening.smact_validity("Al3Li", include_alloys=False)
        self.assertFalse(invalid_result)

    def test_smact_validity_parameters(self):
        """
        Test smact_validity with different parameter combinations
        """
        # Test with pauling test disabled
        result_no_pauling = smact.screening.smact_validity("NaCl", use_pauling_test=False)
        self.assertTrue(result_no_pauling)

        # Test with metallicity check
        result_metallicity = smact.screening.smact_validity("Fe", check_metallicity=True, metallicity_threshold=0.1)
        self.assertTrue(result_metallicity)

        # Test with different oxidation states set
        result_ox_states = smact.screening.smact_validity("NaCl", oxidation_states_set="icsd16")
        self.assertTrue(result_ox_states)

    def test_smact_validity_type_error_forced(self):
        """
        Force the except TypeError block in smact_validity's try/except to be triggered.
        """

        original_pauling_test = smact.screening.pauling_test

        def mock_pauling_test(*args, **kwargs):
            raise TypeError("Mock TypeError triggered")

        try:
            smact.screening.pauling_test = mock_pauling_test
            result = smact.screening.smact_validity("NaCl", use_pauling_test=True)
            self.assertIsInstance(result, bool)
            self.assertTrue(result)
        finally:
            smact.screening.pauling_test = original_pauling_test

    def test_smact_validity_noble_gases(self):
        """
        Test that noble gases with no oxidation states return False cleanly
        (no TypeError) regardless of which oxidation_states_set is used.
        """
        for ox_set in ["smact14", "icsd16", "icsd24", "pymatgen_sp"]:
            self.assertFalse(
                smact.screening.smact_validity("NeF2", oxidation_states_set=ox_set),
                f"NeF2 should be False with oxidation_states_set={ox_set}",
            )
        # Xe and Kr have known oxidation states and should still pass
        self.assertTrue(smact.screening.smact_validity("XeF2"))
        self.assertTrue(smact.screening.smact_validity("KrF2"))

    def test_smact_validity_special_cases(self):
        """
        Test special cases in smact_validity:
        - Single elements
        - Alloys
        - High metallicity
        """
        # Single element should be valid
        self.assertTrue(smact.screening.smact_validity("Fe"))

        # Alloy should be valid when include_alloys=True
        self.assertTrue(smact.screening.smact_validity("FeAl", include_alloys=True))

        # High metallicity check
        self.assertTrue(
            smact.screening.smact_validity(
                "FeNi",
                check_metallicity=True,
                metallicity_threshold=0.5,
                include_alloys=False,  # Test metallicity path specifically
            )
        )

    def test_smact_validity_oxidation_states(self):
        """
        Test various oxidation state sets in smact_validity
        """
        # Test different oxidation state sets
        oxidation_sets = ["smact14", "icsd16", "icsd24", "pymatgen_sp"]
        for ox_set in oxidation_sets:
            self.assertTrue(
                smact.screening.smact_validity("NaCl", oxidation_states_set=ox_set),
                f"Failed with oxidation state set: {ox_set}",
            )

        # Test with file path

        files_dir = join(dirname(realpath(__file__)), "files")
        test_ox_states = join(files_dir, "test_oxidation_states.txt")
        self.assertTrue(smact.screening.smact_validity("NaCl", oxidation_states_set=test_ox_states))

    def test_smact_validity_mixed_valence(self):
        """Test mixed valence handling in smact_validity."""
        # Fe3O4: fails without mixed_valence (no single Fe oxidation state works)
        self.assertFalse(smact.screening.smact_validity("Fe3O4"))
        # Fe3O4: passes with mixed_valence=True (Fe2+ + 2*Fe3+ + 4*O2-)
        self.assertTrue(smact.screening.smact_validity("Fe3O4", mixed_valence=True))

        # Mn3O4 (hausmannite): Mn2+Mn3+2O4 — second real mixed-valence compound
        self.assertTrue(smact.screening.smact_validity("Mn3O4", mixed_valence=True))

        # Compounds with no mixed-valence elements: flag has no effect on valid compounds
        self.assertTrue(smact.screening.smact_validity("NaCl", mixed_valence=True))

        # Compound with no mixed-valence elements that fails regardless (Li only has +1, can't balance)
        self.assertFalse(smact.screening.smact_validity("LiF2", mixed_valence=True))

        # Explicit confirmation of default behaviour (mixed_valence=False)
        self.assertFalse(smact.screening.smact_validity("Fe3O4", mixed_valence=False))

    def test_smact_validity_error_handling(self):
        """
        Test error handling in smact_validity
        """
        # Test with invalid oxidation states set
        with pytest.raises(ValueError):
            smact.screening.smact_validity("NaCl", oxidation_states_set="invalid_set")

        # Test that wiki set gives warning
        with pytest.warns(UserWarning, match=r"This set of oxidation states is from Wikipedia"):
            smact.screening.smact_validity("NaCl", oxidation_states_set="wiki")

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

    # --------- Additional Coverage Tests (for untested lines in screening) ---------

    def test_smact_validity_oxidation_states_file(self):
        """
        Test that providing a valid file path triggers the file-based combos branch.
        """

        self.assertTrue(exists(TEST_OX_STATES), "TEST_OX_STATES must exist for this test.")
        # Test just the validity
        self.assertTrue(smact.screening.smact_validity("NaCl", oxidation_states_set=TEST_OX_STATES))
