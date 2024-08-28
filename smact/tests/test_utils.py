import unittest

from smact.utils.composition import parse_formula


class TestComposition(unittest.TestCase):
    """Test composition utilities"""

    def test_parse_formula(self):
        formulas = ["Li10GeP2S12", "Mg0.5O0.5", "CaMg(CO3)2"]

        LGPS = parse_formula(formulas[0])
        self.assertIsInstance(LGPS, dict)
        for el_sym, ammt in LGPS.items():
            self.assertIsInstance(el_sym, str)
            self.assertIsInstance(ammt, float)
        self.assertEqual(LGPS["Li"], 10)
        self.assertEqual(LGPS["Ge"], 1)
        self.assertEqual(LGPS["P"], 2)
        self.assertEqual(LGPS["S"], 12)

        MgO = parse_formula(formulas[1])
        self.assertIsInstance(MgO, dict)
        self.assertEqual(MgO["Mg"], 0.5)
        self.assertEqual(MgO["O"], 0.5)

        dolomite = parse_formula(formulas[2])
        self.assertIsInstance(dolomite, dict)
        self.assertEqual(dolomite["Ca"], 1)
        self.assertEqual(dolomite["Mg"], 1)
        self.assertEqual(dolomite["C"], 2)
        self.assertEqual(dolomite["O"], 6)
