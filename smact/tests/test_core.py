#!/usr/bin/env python

import os
import unittest

from pymatgen.core.periodic_table import Specie

import smact
import smact.distorter
import smact.lattice
import smact.lattice_parameters
import smact.oxidation_states
import smact.screening
from smact import Species
from smact.builder import wurtzite
from smact.properties import band_gap_Harrison, compound_electroneg


class TestSequenceFunctions(unittest.TestCase):

    # ---------------- TOP-LEVEL ----------------

    def test_Element_class_Pt(self):
        Pt = smact.Element("Pt")
        self.assertEqual(Pt.name, "Platinum")
        self.assertEqual(Pt.ionpot, 8.95883)
        self.assertEqual(Pt.number, 78)
        self.assertEqual(Pt.dipol, 44.00)

    def test_ordered_elements(self):
        self.assertEqual(
            smact.ordered_elements(65, 68), ["Tb", "Dy", "Ho", "Er"]
        )
        self.assertEqual(smact.ordered_elements(52, 52), ["Te"])

    def test_element_dictionary(self):
        newlist = ["O", "Rb", "W"]
        dictionary = smact.element_dictionary(newlist)
        self.assertEqual(dictionary["O"].crustal_abundance, 461000.0)
        self.assertEqual(dictionary["Rb"].oxidation_states, [-1, 1])
        self.assertEqual(dictionary["W"].name, "Tungsten")
        self.assertTrue("Rn" in smact.element_dictionary())

    def test_are_eq(self):
        self.assertTrue(
            smact.are_eq(
                [1.00, 2.00, 3.00], [1.001, 1.999, 3.00], tolerance=1e-2
            )
        )
        self.assertFalse(
            smact.are_eq([1.00, 2.00, 3.00], [1.001, 1.999, 3.00])
        )

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
            compound_electroneg(
                elements=["Cu", "Zn"], stoichs=[0.5, 0.5], source="Pauling"
            ),
            5.0638963259,
        )

    def test_harrison_gap_MgCl(self):
        self.assertAlmostEqual(
            band_gap_Harrison("Mg", "Cl", verbose=False, distance=2.67),
            3.545075110572662,
        )

    # ---------------- BUILDER ----------------

    def test_builder_ZnS(self):
        ZnS, sys_ZnS = wurtzite(["Zn", "S"])
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
        self.assertFalse(
            smact.screening.pauling_test(
                (-2, +2), (Sn.pauling_eneg, S.pauling_eneg)
            )
        )
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
            smact.screening.eneg_states_test(
                [1, 3, -1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg]
            )
        )
        self.assertFalse(
            smact.screening.eneg_states_test(
                [-1, 3, 1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg]
            )
        )

    def test_eneg_states_test_alternate(self):
        Na, Fe, Cl = (smact.Element(label) for label in ("Na", "Fe", "Cl"))
        self.assertTrue(
            smact.screening.eneg_states_test_alternate(
                [1, 3, -1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg]
            )
        )
        self.assertFalse(
            smact.screening.eneg_states_test_alternate(
                [-1, 3, 1], [Na.pauling_eneg, Fe.pauling_eneg, Cl.pauling_eneg]
            )
        )

    def test_eneg_states_test_threshold(self):
        self.assertFalse(
            smact.screening.eneg_states_test_threshold(
                [1, -1], [1.83, 1.82], threshold=0
            )
        )
        self.assertTrue(
            smact.screening.eneg_states_test_threshold(
                [1, -1], [1.83, 1.82], threshold=0.1
            )
        )

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
        self.assertEqual(
            smact.screening.ml_rep_generator(["Pb", "O"], [1, 2]), PbO2_ml
        )
        self.assertEqual(
            smact.screening.ml_rep_generator([Pb, O], [1, 2]), PbO2_ml
        )

    def test_smact_filter(self):
        Na, Fe, Cl = (smact.Element(label) for label in ("Na", "Fe", "Cl"))
        self.assertEqual(
            smact.screening.smact_filter([Na, Fe, Cl], threshold=2),
            [
                (("Na", "Fe", "Cl"), (1, -1, -1), (2, 1, 1)),
                (("Na", "Fe", "Cl"), (1, 1, -1), (1, 1, 2)),
            ],
        )
        self.assertEqual(
            len(smact.screening.smact_filter([Na, Fe, Cl], threshold=8)), 77
        )

    # ---------------- Lattice ----------------
    def test_Lattice_class(self):
        site_A = smact.lattice.Site([0, 0, 0], -1)
        site_B = smact.lattice.Site([0.5, 0.5, 0.5], [+2, +3])
        test_lattice = smact.lattice.Lattice([site_A, site_B], space_group=221)
        self.assertEqual(test_lattice.sites[0].oxidation_states, -1)
        self.assertEqual(test_lattice.sites[1].position, [0.5, 0.5, 0.5])

    # ---------- Lattice parameters -----------
    def test_lattice_parameters(self):
        perovskite = smact.lattice_parameters.cubic_perovskite(
            [1.81, 1.33, 1.82]
        )
        wurtz = smact.lattice_parameters.wurtzite([1.81, 1.33])
        self.assertAlmostEqual(perovskite[0], 6.3)
        self.assertAlmostEqual(perovskite[1], 6.3)
        self.assertAlmostEqual(perovskite[3], 90)
        self.assertAlmostEqual(wurtz[0], 5.13076)
        self.assertAlmostEqual(wurtz[2], 8.3838)

    # ---------- Lattice parameters -----------
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
