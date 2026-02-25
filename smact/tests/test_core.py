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
            0.0,  # element 103 (Lr)
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
        # Also verify the explicit oxidation_states_set path (the reported bug case)
        self.assertTrue(smact.screening.smact_validity("XeF2", oxidation_states_set="icsd24"))
        self.assertTrue(smact.screening.smact_validity("KrF2", oxidation_states_set="icsd24"))

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

    # --------- Properties branch coverage ---------

    def test_eneg_mulliken_branches(self):
        """eneg_mulliken: str input, Element input, and invalid type."""
        from smact.properties import eneg_mulliken

        # str input (line 32)
        val_str = eneg_mulliken("Fe")
        self.assertIsInstance(val_str, float)

        # Element input (line 36)
        val_el = eneg_mulliken(smact.Element("Fe"))
        self.assertAlmostEqual(val_str, val_el)

        # Invalid type raises TypeError (lines 33-34)
        with pytest.raises(TypeError):
            eneg_mulliken(42)

    def test_harrison_gap_verbose(self):
        """band_gap_Harrison with verbose=True hits lines 89-92."""
        import io
        import sys

        buf = io.StringIO()
        sys.stdout = buf
        try:
            band_gap_Harrison("Mg", "Cl", verbose=True, distance=2.67)
        finally:
            sys.stdout = sys.__stdout__
        out = buf.getvalue()
        self.assertIn("V1_bar", out)
        self.assertIn("V2", out)
        self.assertIn("alpha_m", out)
        self.assertIn("V3", out)

    def test_compound_electroneg_element_objects(self):
        """compound_electroneg with Element objects (lines 133-134) and verbose (155, 167)."""
        Fe = smact.Element("Fe")
        O = smact.Element("O")

        # Element list input (lines 133-134) + Mulliken (line 145)
        val = compound_electroneg(elements=[Fe, O], stoichs=[1, 1], source="Mulliken")
        self.assertIsInstance(val, float)

        # Invalid element type raises TypeError (lines 135-136)
        with pytest.raises(TypeError):
            compound_electroneg(elements=[42, 43], stoichs=[1, 1])

        # Invalid source raises ValueError (line 150)
        with pytest.raises(ValueError):
            compound_electroneg(elements=["Fe", "O"], stoichs=[1, 1], source="BadSource")

        # verbose=True covers lines 155 and 167
        import io
        import sys

        buf = io.StringIO()
        sys.stdout = buf
        try:
            compound_electroneg(elements=["Fe", "O"], stoichs=[1, 1], verbose=True)
        finally:
            sys.stdout = sys.__stdout__
        out = buf.getvalue()
        self.assertIn("Electronegativities", out)
        self.assertIn("Geometric mean", out)

    # --------- Builder branch coverage ---------

    def test_cubic_perovskite(self):
        """cubic_perovskite was never tested; covers builder.py lines 41-58."""
        from smact.builder import cubic_perovskite

        lattice, system = cubic_perovskite(["Ba", "Ti", "O"])
        self.assertEqual(len(lattice.sites), 5)  # 1 Ba + 1 Ti + 3 O
        self.assertIsNotNone(system)

    # --------- Oxidation-states branch coverage ---------

    def test_oxidation_states_error_branches(self):
        """OxidationStateProbabilityFinder error paths (lines 80, 88, 151, 161, 182)."""
        ox = smact.oxidation_states.OxidationStateProbabilityFinder()

        # Both positive → ValueError (line 80)
        with pytest.raises(ValueError, match="cation and one anion"):
            ox._generate_lookup_key(Species("Fe", 3), Species("Cu", 2))

        # Both negative → ValueError (line 80)
        with pytest.raises(ValueError, match="cation and one anion"):
            ox._generate_lookup_key(Species("O", -2), Species("S", -2))

        # Species not in table → NameError (line 88)
        # Fe+5 is a valid smact Species but is absent from the Hautier probability table
        with pytest.raises(NameError):
            ox._generate_lookup_key(Species("Fe", 5), Species("O", -2))

        # Invalid input type → TypeError (line 151)
        with pytest.raises(TypeError):
            ox.compound_probability("NaCl_string")

        # No cations in list → ValueError (line 161)
        with pytest.raises(ValueError, match="No cations"):
            ox.compound_probability([Species("O", -2), Species("S", -2)])

        # Unknown module attribute → AttributeError (line 182)
        with pytest.raises(AttributeError):
            _ = smact.oxidation_states.does_not_exist

    def test_oxidation_states_pymatgen_species_input(self):
        """compound_probability with pymatgen Species list (lines 134-138, 142)."""
        from pymatgen.core import Lattice as PmgLattice
        from pymatgen.core import Structure
        from pymatgen.core.periodic_table import Specie as PmgSpecie

        ox = smact.oxidation_states.OxidationStateProbabilityFinder()

        # List of pymatgen Species (lines 134-136)
        result = ox.compound_probability([PmgSpecie("Fe", 3), PmgSpecie("O", -2)])
        self.assertIsInstance(result, float)

        # List with mixed/invalid pymatgen species (line 138) — a list mixing types
        with pytest.raises(TypeError):
            ox.compound_probability([PmgSpecie("Fe", 3), "not_a_species"])

        # pymatgen Structure without oxidation states → TypeError (line 142)
        struct = Structure.from_spacegroup(
            "Fm-3m", PmgLattice.cubic(5.6), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        )
        with pytest.raises(TypeError, match="oxidation states"):
            ox.compound_probability(struct)

    # --------- Screening branch coverage ---------

    def test_pauling_test_threshold_branches(self):
        """pauling_test lines 139 (repeat+threshold) and 150 (no-repeat+threshold)."""
        Sn, S = smact.Element("Sn"), smact.Element("S")

        # Both repeats=True with threshold != 0 → line 139
        result = smact.screening.pauling_test(
            (+2, -2),
            (Sn.pauling_eneg, S.pauling_eneg),
            repeat_anions=True,
            repeat_cations=True,
            threshold=0.5,
        )
        self.assertIsInstance(result, bool)

        # repeat_anions=False, repeat_cations=False, distinct symbols, threshold != 0
        # → _no_repeats line 179 returns True → pauling_test line 150
        result2 = smact.screening.pauling_test(
            (-2, +2),
            (S.pauling_eneg, Sn.pauling_eneg),
            symbols=("S", "Sn"),
            repeat_anions=False,
            repeat_cations=False,
            threshold=0.5,
        )
        self.assertIsInstance(result2, bool)

    def test_no_repeats_line_179(self):
        """_no_repeats line 179: repeat_anions=False AND repeat_cations=False → len check."""
        Sn, S = smact.Element("Sn"), smact.Element("S")

        # Both repeat flags False, distinct symbols → _no_repeats line 179 → True → line 148
        result = smact.screening.pauling_test(
            (-2, +2),
            (S.pauling_eneg, Sn.pauling_eneg),
            symbols=("S", "Sn"),
            repeat_anions=False,
            repeat_cations=False,
            threshold=0.0,
        )
        self.assertIsInstance(result, bool)

    def test_pauling_test_old_edge_cases(self):
        """pauling_test_old: None in paul (236), all-positive (263), equal max/min (265)."""
        Sn, S = smact.Element("Sn"), smact.Element("S")

        # None in pauling array → returns False immediately (line 236)
        with pytest.warns(DeprecationWarning):
            result = smact.screening.pauling_test_old(
                (+2, -2), (None, S.pauling_eneg), symbols=("Sn", "S")
            )
        self.assertFalse(result)

        # All-positive oxidation states → len(negative)==0 → False (line 263)
        with pytest.warns(DeprecationWarning):
            result2 = smact.screening.pauling_test_old(
                (+2, +3), (Sn.pauling_eneg, 2.0), symbols=("Sn", "Fe")
            )
        self.assertFalse(result2)

        # max(positive) == min(negative) → False (line 265)
        with pytest.warns(DeprecationWarning):
            result3 = smact.screening.pauling_test_old(
                (+2, -2), (1.8, 1.8), symbols=("Sn", "S")
            )
        self.assertFalse(result3)

    def test_ml_rep_generator_no_stoichs(self):
        """ml_rep_generator without explicit stoichs triggers default-stoichs path (line 403)."""
        result = smact.screening.ml_rep_generator(["Na", "Cl"])
        self.assertEqual(len(result), 103)
        self.assertAlmostEqual(sum(result), 1.0)
        # Each element contributes equally when stoichs default to 1
        Na_idx = int(smact.Element("Na").number) - 1
        Cl_idx = int(smact.Element("Cl").number) - 1
        self.assertAlmostEqual(result[Na_idx], 0.5)
        self.assertAlmostEqual(result[Cl_idx], 0.5)

    def test_smact_filter_wiki_warning(self):
        """smact_filter with wiki oxidation_states_set emits a warning (line 487)."""
        Na, Cl = smact.Element("Na"), smact.Element("Cl")
        with pytest.warns(UserWarning, match="Wikipedia"):
            smact.screening.smact_filter([Na, Cl], threshold=2, oxidation_states_set="wiki")

    def test_smact_filter_invalid_set_raises(self):
        """smact_filter with unknown set string raises ValueError (line 494)."""
        Na, Cl = smact.Element("Na"), smact.Element("Cl")
        with pytest.raises(ValueError):
            smact.screening.smact_filter([Na, Cl], oxidation_states_set="nonexistent_set")

    def test_smact_filter_missing_oxidation_states_raises(self):
        """smact_filter raises ValueError when element has no ox states in chosen set (line 501)."""
        # Argon has no oxidation states in icsd24
        Ar, Na = smact.Element("Ar"), smact.Element("Na")
        with pytest.raises(ValueError, match="No oxidation states"):
            smact.screening.smact_filter([Na, Ar], threshold=2, oxidation_states_set="icsd24")

    def test_smact_validity_filter_params_warning(self):
        """smact_validity warns when filtering params used with explicit oxidation_states_set (line 561)."""
        with pytest.warns(UserWarning, match="include_zero"):
            smact.screening.smact_validity("NaCl", oxidation_states_set="icsd24", include_zero=True)

    def test_smact_validity_element_not_in_consensus_set(self):
        """smact_validity returns False when element absent from consensus ICSD24 filter (line 609)."""
        # Noble gas with no ICSD entries: Ar not expected in consensus filter
        result = smact.screening.smact_validity("ArNa", oxidation_states_set=None)
        self.assertFalse(result)

    # --------- distorter coverage ---------

    def test_distorter_functions(self):
        """distorter get_sg, build_sub_lattice, get_inequivalent_sites, make_substitution (lines 44-47, 69-87, 105-114, 132-140).

        Newer spglib dropped direct ASE Atoms support (returns None), so we patch
        spglib.get_spacegroup to return the canonical Fm-3m string for NaCl.
        """
        from unittest.mock import patch

        from ase.spacegroup import Spacegroup
        from ase.spacegroup import crystal as ase_crystal

        import smact.distorter as distorter

        # NaCl rocksalt: spacegroup 225 (Fm-3m)
        nacl = ase_crystal(
            ["Na", "Cl"],
            [(0, 0, 0), (0.5, 0.5, 0.5)],
            spacegroup=225,
            cellpar=[5.6, 5.6, 5.6, 90, 90, 90],
        )

        with patch.object(distorter.spglib, "get_spacegroup", return_value="Fm-3m (225)"):
            # get_sg (lines 44-47)
            sg = distorter.get_sg(nacl)
            self.assertIsInstance(sg, Spacegroup)
            self.assertEqual(sg.no, 225)

            # build_sub_lattice (lines 132-140) — does not use spglib
            sub_lat = distorter.build_sub_lattice(nacl, "Na")
            self.assertGreater(len(sub_lat), 0)
            self.assertEqual(len(sub_lat[0]), 3)

            # get_inequivalent_sites (lines 69-87) — calls get_sg internally
            ineq = distorter.get_inequivalent_sites(sub_lat, nacl)
            self.assertGreater(len(ineq), 0)

            # make_substitution (lines 105-114) — does not use spglib
            new_lattice = distorter.make_substitution(nacl, sub_lat[0], "K")
            self.assertIn("K", new_lattice.get_chemical_symbols())

            # get_inequivalent_sites line 77: duplicate site → new_site = False
            # Passing the same site twice makes the second pass find a match
            ineq_dup = distorter.get_inequivalent_sites([sub_lat[0], sub_lat[0]], nacl)
            self.assertEqual(len(ineq_dup), 1)

    def test_smact_validity_mixed_valence_blowup_warning(self):
        """smact_validity warns and returns False when MV expansion exceeds 1M combinations (lines 645-650).

        Fe8O17: gcd(8,17)=1 so the formula is irreducible → stoich=(8,17).
        8*x + 17*(-2) = 0 requires x=4.25 (non-integer), so standard check fails.
        projected = len(Fe_smact14_ox)^8 * len(O_smact14_ox) = 8^8 * 2 ≈ 33M > 1M.
        """
        with pytest.warns(UserWarning, match="too many combinations"):
            result = smact.screening.smact_validity(
                "Fe8O17", mixed_valence=True, oxidation_states_set="smact14"
            )
        self.assertFalse(result)
