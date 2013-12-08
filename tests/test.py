#!/usr/bin/env python

import unittest
from compound_electroneg import compound_electroneg

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass
    
    def test_compound_eneg_brass(self):
        self.assertEqual(compound_electroneg(
            elements="Cu Zn", stoichs="0.5 0.5"),
            4.5878238674779128)

if __name__ == '__main__':
    unittest.main()
