import unittest

import smact
from smact.dopant_prediction import doper
from smact.structure_prediction import mutation, utilities


class DopantTest(unittest.TestCase):
    def test_dopant_prediction(self):
        test = doper.Doper(("Cd2+", "O2-"))

        test_cation = "Cd2+"
        cation, charge = utilities.parse_spec(test_cation)

        # Assert: Length of the list and return type (dictionary: list)
        self.assertIs(type(test.get_dopants()), dict)
        for ls in test.get_dopants().values():
            self.assertEqual(len(test), len(ls))
            self.assertIs(type(ls), list)

        # Assert: (cation) higher charges for n-type and lower charges for p-type
        n_sub_list = test.get_dopants().get("n-type cation substitutions")
        p_sub_list = test.get_dopants().get("p-type cation substitutions")
        for n_atom, p_atom in zip(n_sub_list, p_sub_list):
            self.assertGreater(utilities.parse_spec(n_atom[0])[1], charge)
            self.assertLess(utilities.parse_spec(p_atom[0])[1], charge)

if __name__ == "__main__":
    TestLoader = unittest.TestLoader()
    DoperTests = unittest.TestSuite()
    DoperTests.addTests(TestLoader.loadTestsFromTestCase(DopantTest))

    runner = unittest.TextTestRunner()
    result = runner.run(DoperTests)
