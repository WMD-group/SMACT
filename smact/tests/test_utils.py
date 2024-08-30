import unittest

from pymatgen.core import Composition

from smact import Element
from smact.screening import smact_filter
from smact.utils.composition import comp_maker, formula_maker, parse_formula


class TestComposition(unittest.TestCase):
    """Test composition utilities"""

    def setUp(self) -> None:
        self.mock_filter_output = [
            (("Fe", "O"), (2, -2), (1, 1)),
            (("Fe", "O"), (1, 1)),
            (("Fe", "Fe", "O"), (2, 3, -2), (1, 2, 4)),
        ]
        self.smact_filter_output = smact_filter(
            els=[Element("Li"), Element("Ge"), Element("P"), Element("S")],
            stoichs=[[10], [1], [2], [12]],
        )

    def test_parse_formula(self):
        """Test the parse_formula function"""

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

    def test_comp_maker(self):
        """Test the comp_maker function"""
        comp1 = comp_maker(self.mock_filter_output[0])
        comp2 = comp_maker(self.mock_filter_output[1])
        comp3 = comp_maker(self.mock_filter_output[2])
        comp4 = comp_maker(self.smact_filter_output[1])
        for comp in [comp1, comp2, comp3, comp4]:
            self.assertIsInstance(comp, Composition)
        self.assertEqual(Composition("FeO"), comp2)
        self.assertEqual(Composition({"Fe2+": 1, "O2-": 1}), comp1)
        self.assertEqual(Composition({"Fe2+": 1, "Fe3+": 2, "O2-": 4}), comp3)
        self.assertEqual(
            Composition({"Li+": 10, "Ge4+": 1, "P5+": 2, "S2-": 12}), comp4
        )

    def test_formula_maker(self):
        """Test the formula_maker function"""
        form1 = formula_maker(self.mock_filter_output[0])
        form2 = formula_maker(self.mock_filter_output[1])
        form3 = formula_maker(self.mock_filter_output[2])
        form4 = formula_maker(self.smact_filter_output[1])
        self.assertEqual(form1, "FeO")
        self.assertEqual(form2, "FeO")
        self.assertEqual(form1, form2)
        self.assertEqual(form3, "Fe3O4")
        self.assertEqual(form4, "Li10Ge(PS6)2")
