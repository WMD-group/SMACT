#!/usr/bin/env python

import unittest
from smact.properties.compound_electroneg import compound_electroneg
from smact.builder import wurtzite
from smact.lattice import *
import os
import smact.core

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_Element_class_Pt(self):
        Pt = smact.core.Element('Pt')
        self.assertEqual(Pt.name,'Platinum')
        self.assertEqual(Pt.ionpot,8.9588)
        self.assertEqual(Pt.number,78)

#TODO(AJJ): Implement warnings in a testable way. Have a Pu testcase.
    def test_compound_eneg_brass(self):
        self.assertEqual(compound_electroneg(
            elements=["Cu","Zn"], stoichs=[0.5, 0.5]),
            4.5878238674779128)

    def test_builder_ZnS(self):
	ZnS = wurtzite(['Zn','S'])
	self.assertEqual((ZnS.sites[0].position[2]),0)
	self.assertEqual((ZnS.sites[0].position[0]),2./3.)
'''    

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
        f = open(smact_directory + 'data/oxidationstates.data','r')
        lines = f.readlines()
        f.close()
        for line in lines:
            inp = line.split()
            for item in search_space:
                if inp[0].startswith(item):
                    key = inp[0]
                    elements[key] = inp[1]

# Generate list of compositions which satisfy charge neutrality
        perovskite = Lattice(["A","B","C"],[1,1,3],[[1,2],[2,3,4],[-1,-2]])
        perovskite_compositions = possible_compositions(perovskite, elements)
        self.assertEqual(perovskite_compositions[3][2],'C-2')
        self.assertEqual(perovskite_compositions[12][2],'S-2')
        self.assertEqual(perovskite_compositions[20][1],'Ge4')
        self.assertEqual(perovskite_compositions[32][0],'S2')
#--------------------------------------------------------------------------------------------------------
'''

if __name__ == '__main__':
    unittest.main()
