#!/usr/bin/env python

import unittest
from smact.properties.compound_electroneg import compound_electroneg
import smact.smact_core as smact_core
from smact.smact_builder import wurtzite

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_Element_class_Pt(self):
        Pt = smact_core.Element('Pt')
        self.assertEqual(Pt.name,'Platinum')
        self.assertEqual(Pt.ionpot,8.9588)

#TODO(AJJ): Implement warnings in a testable way. Have a Pu testcase.
    
    def test_compound_eneg_brass(self):
        self.assertEqual(compound_electroneg(
            elements=["Cu","Zn"], stoichs=[0.5, 0.5]),
            4.5878238674779128)

    def test_builder_ZnS(self):
	ZnS = wurtzite(['Zn','S'])
	self.assertEqual(round(ZnS.positions[0,0],2),1.5)
	self.assertEqual(round(ZnS.positions[1,1],2),1.73)
	self.assertEqual(round(ZnS.positions[2,2],2),3.75)
	self.assertEqual(round(ZnS.cell[0,0],2),3.0)
	self.assertEqual(round(ZnS.cell[1,0],2),-1.5)

if __name__ == '__main__':
    unittest.main()
