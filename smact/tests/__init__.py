import unittest

import smact.tests.test
from smact.tests.test_structure import CationMutatorTest, StructureDBTest, StructureTest

if __name__ == "__main__":
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromModule(smact.tests.test))
    suite.addTests(loader.loadTestsFromTestCase(StructureTest))
    suite.addTests(loader.loadTestsFromTestCase(StructureDBTest))
    suite.addTests(loader.loadTestsFromTestCase(CationMutatorTest))
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
