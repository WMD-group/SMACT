#!/usr/bin/env python

import unittest
from smact.properties import compound_electroneg
from smact.builder import wurtzite
import smact.lattice
import smact

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass

    ################ TOP-LEVEL ################
    
    def test_Element_class_Pt(self):
        Pt = smact.Element('Pt')
        self.assertEqual(Pt.name,'Platinum')
        self.assertEqual(Pt.ionpot,8.9588)
        self.assertEqual(Pt.number,78)

    def test_ordered_elements(self):
        self.assertEqual(
            smact.ordered_elements(65, 68),
            ['Tb', 'Dy', 'Ho', 'Er']
            )
        self.assertEqual(
            smact.ordered_elements(52,52),
            ['Te']
            )
    def test_element_dictionary(self):
        newlist = ['O','Rb','W']
        dictionary = smact.element_dictionary(newlist)
        self.assertEqual(dictionary['O'].crustal_abundance, 461000.0)
        self.assertEqual(dictionary['Rb'].oxidation_states,[-1, 1])
        self.assertEqual(dictionary['W'].name,'Tungsten')
        self.assertTrue('Rn' in smact.element_dictionary())

    def test_are_eq(self):
        self.assertTrue(
            smact.are_eq([1.00, 2.00, 3.00],
                         [1.001, 1.999, 3.00],
                         tolerance=1e-2)
            )
        self.assertFalse(
            smact.are_eq([1.00, 2.00, 3.00],
                         [1.001, 1.999, 3.00])
            )

    def test_gcd_recursive(self):
        self.assertEqual(
            smact._gcd_recursive(
                4, 12, 10, 32
                ),
            2)
        self.assertEqual(
            smact._gcd_recursive(
                15, 4
                ),
            1)

    def test_isneutral(self):
        self.assertTrue(
            smact._isneutral(
            (-2,1), (2, 4)
            ))
        self.assertFalse(
            smact._isneutral(
            (4,-1), (1, 3)
            ))
        
    def test_charge_neutrality_ternary(self):
        ox = [1,-2,1]
        is_neutral, neutral_combos = smact.charge_neutrality(ox)
        self.assertTrue((is_neutral))
        self.assertEqual(len(neutral_combos),9)
        self.assertTrue((3, 2, 1) in neutral_combos)

    def test_pauling_test(self):
        Sn, S = (smact.Element(label) for label in ('Sn', 'S'))
        self.assertTrue(smact.pauling_test(
            (+2, -2), (Sn.pauling_eneg, S.pauling_eneg)
            ))
        self.assertFalse(smact.pauling_test(
            (-2, +2), (Sn.pauling_eneg, S.pauling_eneg)
            ))
                        
    ################ PROPERTIES ################

    def test_compound_eneg_brass(self):
        self.assertEqual(compound_electroneg(
            elements=["Cu","Zn"], stoichs=[0.5, 0.5]),
            4.5878238674779128)

    ################ BUILDER ################
        
    def test_builder_ZnS(self):
        ZnS, sys_ZnS = wurtzite(['Zn','S'])
        self.assertEqual((ZnS.sites[0].position[2]),0)
        self.assertEqual((ZnS.sites[0].position[0]),2./3.)



''''
    def test_compound_library(self):
        """A compound test, this covers builder, oxidation_data and the possible_compositions functions"""
        this_directory = os.path.dirname(__file__)
        if this_directory:
            this_directory = this_directory + '/'
        smact_directory = this_directory + '../smact/'

# Generate a dictionary elements, form the dataset oxidationstates.data
# Dictionary contains elements and their oxidation states
# Reduce the regions of the periodic table to visit, by using search_space
        search_space = {'Ti','O-','Sn','F-','C','Sr','Mg','Cu','Li','S','Si','Ge'}
# Get the list of possible constituant elements
        elements = {}
        f = open(path.join(smact.data_directory, 'oxidationstates.data'), 'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            inp = line.split()
            for item in search_space:
                if inp[0].startswith(item):
                    key = inp[0]
                    elements[key] = inp[1]

# Generate list of compositions which satisfy charge neutrality
        perovskite = smact.lattice.Lattice(["A","B","C"],[1,1,3],[[1,2],[2,3,4],[-1,-2]])
        perovskite_compositions = smact.lattice.possible_compositions(perovskite, elements)
        self.assertEqual(perovskite_compositions[3][2],'C-2')
        self.assertEqual(perovskite_compositions[12][2],'S-2')
        self.assertEqual(perovskite_compositions[20][1],'Ge4')
        self.assertEqual(perovskite_compositions[32][0],'S2')
#--------------------------------------------------------------------------------------------------------
'''

if __name__ == '__main__':
    unittest.main()
