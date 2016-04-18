#!/usr/bin/env python

import unittest
import smact
from smact.properties import compound_electroneg
from smact.builder import wurtzite
import smact.screening
import smact.lattice


class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass

    # ---------------- TOP-LEVEL ----------------

    def test_Element_class_Pt(self):
        Pt = smact.Element('Pt')
        self.assertEqual(Pt.name, 'Platinum')
        self.assertEqual(Pt.ionpot, 8.95883)
        self.assertEqual(Pt.number, 78)

    def test_ordered_elements(self):
        self.assertEqual(
            smact.ordered_elements(65, 68),
            ['Tb', 'Dy', 'Ho', 'Er'])
        self.assertEqual(
            smact.ordered_elements(52, 52),
            ['Te'])

    def test_element_dictionary(self):
        newlist = ['O', 'Rb', 'W']
        dictionary = smact.element_dictionary(newlist)
        self.assertEqual(dictionary['O'].crustal_abundance, 461000.0)
        self.assertEqual(dictionary['Rb'].oxidation_states, [-1, 1])
        self.assertEqual(dictionary['W'].name, 'Tungsten')
        self.assertTrue('Rn' in smact.element_dictionary())

    def test_are_eq(self):
        self.assertTrue(
            smact.are_eq([1.00, 2.00, 3.00],
                         [1.001, 1.999, 3.00],
                         tolerance=1e-2))
        self.assertFalse(
            smact.are_eq([1.00, 2.00, 3.00],
                         [1.001, 1.999, 3.00]))

    def test_gcd_recursive(self):
        self.assertEqual(
            smact._gcd_recursive(4, 12, 10, 32),
            2)
        self.assertEqual(
            smact._gcd_recursive(15, 4),
            1)

    def test_isneutral(self):
        self.assertTrue(
            smact._isneutral((-2, 1), (2, 4)))
        self.assertFalse(
            smact._isneutral((4, -1), (1, 3)))

    def test_neutral_ratios(self):
        ox = [1, -2, 1]
        is_neutral, neutral_combos = smact.neutral_ratios(ox)
        self.assertTrue((is_neutral))
        self.assertEqual(len(neutral_combos), 9)
        self.assertTrue((3, 2, 1) in neutral_combos)

    def test_pauling_test(self):
        Sn, S = (smact.Element(label) for label in ('Sn', 'S'))
        self.assertTrue(smact.screening.pauling_test(
            (+2, -2), (Sn.pauling_eneg, S.pauling_eneg),
            ))
        self.assertFalse(smact.screening.pauling_test(
            (-2, +2), (Sn.pauling_eneg, S.pauling_eneg)
            ))
        self.assertFalse(smact.screening.pauling_test(
            (-2, -2, +2), (S.pauling_eneg, S.pauling_eneg, Sn.pauling_eneg),
            symbols=('S', 'S', 'Sn'), repeat_anions=False
            ))
        self.assertTrue(smact.screening.pauling_test(
            (-2, -2, +2), (S.pauling_eneg, S.pauling_eneg, Sn.pauling_eneg),
            symbols=('S', 'S', 'Sn'), repeat_cations=False
            ))
        self.assertFalse(smact.screening.pauling_test(
            (-2, +2, +2), (S.pauling_eneg, Sn.pauling_eneg, Sn.pauling_eneg),
            symbols=('S', 'Sn', 'Sn'), repeat_cations=False
            ))
        self.assertTrue(smact.screening.pauling_test(
            (-2, +2, +2), (S.pauling_eneg, Sn.pauling_eneg, Sn.pauling_eneg),
            symbols=('S', 'Sn', 'Sn'), repeat_anions=False
            ))

    # ---------------- Properties ----------------

    def test_compound_eneg_brass(self):
        self.assertAlmostEqual(compound_electroneg(
            elements=["Cu", "Zn"], stoichs=[0.5, 0.5],
            source='Pauling'),
            5.0638963259)

    # ---------------- BUILDER ----------------

    def test_builder_ZnS(self):
        ZnS, sys_ZnS = wurtzite(['Zn', 'S'])
        self.assertEqual((ZnS.sites[0].position[2]), 0)
        self.assertEqual((ZnS.sites[0].position[0]), 2./3.)


if __name__ == '__main__':
    unittest.main()
