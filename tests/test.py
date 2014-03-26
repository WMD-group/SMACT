#!/usr/bin/env python

import unittest
from compound_electroneg import compound_electroneg
import smact_core

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

if __name__ == '__main__':
    unittest.main()
