import unittest

import smact.tests.test
from smact.tests.test_structure import CationMutatorTest, StructureDBTest, StructureTest, PredictorTest

if __name__ == "__main__":
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromModule(smact.tests.test))
    suite.addTests(loader.loadTestsFromTestCase(StructureTest))
    suite.addTests(loader.loadTestsFromTestCase(StructureDBTest))
    suite.addTests(loader.loadTestsFromTestCase(CationMutatorTest))
    suite.addTests(loader.loadTestsFromTestCase(PredictorTest))
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
