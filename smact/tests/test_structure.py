"""Test structure prediction module"""

import itertools
import json
import logging
import os
import pickle
import unittest
from contextlib import contextmanager

import numpy as np
import pandas as pd
import pymatgen

import smact
from smact import Species
from smact.structure_prediction.database import StructureDB
from smact.structure_prediction.mutation import CationMutator
from smact.structure_prediction.structure import SmactStructure

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_STRUCT = os.path.join(files_dir, "test_struct")
TEST_POSCAR = os.path.join(files_dir, "test_poscar.txt")
TEST_PY_STRUCT = os.path.join(files_dir, "pymatgen_structure.json")
TEST_LAMBDA_TAB = os.path.join(files_dir, "test_lambda_tab.json")


def generate_test_structure(comp: str) -> bool:
    """Generate a pickled test structure for comparison."""
    poscar_file = os.path.join(files_dir, f"{comp}.txt")
    s = SmactStructure.from_file(poscar_file)

    with open(TEST_STRUCT, 'wb') as f:
        pickle.dump(s, f)

    with open(TEST_POSCAR, 'w') as f:
        f.write(s.as_poscar())

    return True


@contextmanager
def ignore_warnings(logger: logging.Logger) -> int:
    """Ignore logging warnings."""
    log_lvl_buff = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield log_lvl_buff
    logger.setLevel(log_lvl_buff)


class StructureTest(unittest.TestCase):
    """`SmactStructure` testing."""

    TEST_SPECIES = {
      "CaTiO3": [('Ca', 2, 1), ('Ti', 4, 1), ('O', -2, 3)],
      "NaCl": [('Na', 1, 1), ('Cl', -1, 1)],
      "Fe": [('Fe', 0, 1)],
    }

    def assertStructAlmostEqual(self, s1: SmactStructure, s2: SmactStructure, places: int = 7):
        """Assert that two SmactStructures are almost equal.

        Almost equality dependent on how many decimal places the site coordinates
        are equal to.
        """
        must_equal = ["species", "lattice_param"]
        for cond in must_equal:
            self.assertEqual(getattr(s1, cond), getattr(s2, cond))

        self.assertTrue(np.array_equal(s1.lattice_mat, s2.lattice_mat))
        self.assertEqual(list(s1.sites.keys()), list(s2.sites.keys()))

        for si, sj in zip(s1.sites, s2.sites):
            for c1, c2 in zip(si, sj):
                self.assertAlmostEqual(c1, c2, places=places)

    def test_as_poscar(self):
        """Test POSCAR generation."""
        for comp in self.TEST_SPECIES.keys():
            with self.subTest(comp=comp):
                comp_file = os.path.join(files_dir, f"{comp}.txt")
                with open(comp_file, "r") as f:
                    struct = SmactStructure.from_file(comp_file)
                    self.assertEqual(struct.as_poscar(), f.read())

    @staticmethod
    def _gen_empty_structure(species):
        """Generate an empty set of arguments for `SmactStructure` testing."""
        lattice_mat = np.array([[0] * 3] * 3)

        if isinstance(species[0][0], str):
            species_strs = [
              "{ele}{charge}{sign}".format(
                ele=spec[0],
                charge=abs(spec[1]),
                sign='+' if spec[1] >= 0 else '-', ) for spec in species
            ]
        else:
            species_strs = [
              "{ele}{charge}{sign}".format(
                ele=spec[0].symbol,
                charge=abs(spec[0].oxidation),
                sign='+' if spec[0].oxidation >= 0 else '-',
              ) for spec in species
            ]

        sites = {spec: [[]] for spec in species_strs}
        return species, lattice_mat, sites

    def test_from_py_struct(self):
        """Test generation of SmactStructure from a pymatgen Structure."""
        with open(TEST_PY_STRUCT, 'r') as f:
            d = json.load(f)
            py_structure = pymatgen.Structure.from_dict(d)

        with ignore_warnings(smact.structure_prediction.logger):
            s1 = SmactStructure.from_py_struct(py_structure)

        s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))

        self.assertStructAlmostEqual(s1, s2)

    def test_smactStruc_comp_key(self):
        """Test generation of a composition key for `SmactStructure`s."""
        s1 = SmactStructure(*self._gen_empty_structure([('Ba', 2, 2), ('O', -2, 1), ('F', -1, 2)]))
        s2 = SmactStructure(*self._gen_empty_structure([('Fe', 2, 1), ('Fe', 3, 2), ('O', -2, 4)]))

        Ba = Species('Ba', 2)
        O = Species('O', -2)
        F = Species('F', -1)
        Fe2 = Species('Fe', 2)
        Fe3 = Species('Fe', 3)

        s3 = SmactStructure(*self._gen_empty_structure([(Ba, 2), (O, 1), (F, 2)]))
        s4 = SmactStructure(*self._gen_empty_structure([(Fe2, 1), (Fe3, 2), (O, 4)]))

        Ba_2OF_2 = "Ba_2_2+F_2_1-O_1_2-"
        Fe_3O_4 = "Fe_2_3+Fe_1_2+O_4_2-"
        self.assertEqual(s1.composition(), Ba_2OF_2)
        self.assertEqual(s2.composition(), Fe_3O_4)
        self.assertEqual(s3.composition(), Ba_2OF_2)
        self.assertEqual(s4.composition(), Fe_3O_4)

    def test_smactStruc_from_file(self):
        """Test the `from_file` method of `SmactStructure`."""
        with open(TEST_STRUCT, 'rb') as f:
            s1 = pickle.load(f)

        s2 = SmactStructure.from_file(TEST_POSCAR)

        self.assertEqual(s1, s2)

    def test_equality(self):
        """Test equality determination of `SmactStructure`."""
        struct_files = [os.path.join(files_dir, f"{x}.txt") for x in ["CaTiO3", "NaCl"]]
        CaTiO3 = SmactStructure.from_file(struct_files[0])
        NaCl = SmactStructure.from_file(struct_files[1])

        with self.subTest(msg="Testing equality of same object."):
            self.assertEqual(CaTiO3, CaTiO3)

        with self.subTest(msg="Testing inequality of different types."):
            self.assertNotEqual(CaTiO3, "CaTiO3")

        with self.subTest(msg="Testing inequality of different objects."):
            self.assertNotEqual(CaTiO3, NaCl)

    def test_ele_stoics(self):
        """Test acquiring element stoichiometries."""
        s1 = SmactStructure(*self._gen_empty_structure([('Fe', 2, 1), ('Fe', 3, 2), ('O', -2, 4)]))
        s1_stoics = {'Fe': 3, 'O': 4}
        s2 = SmactStructure(*self._gen_empty_structure([('Ba', 2, 2), ('O', -2, 1), ('F', -1, 2)]))
        s2_stoics = {'Ba': 2, 'O': 1, 'F': 2}

        for test, expected in [(s1, s1_stoics), (s2, s2_stoics)]:
            with self.subTest(species=test.species):
                self.assertEqual(SmactStructure._get_ele_stoics(test.species), expected)

    @unittest.skipUnless(os.environ.get("MPI_KEY"), "requires MPI key to be set.")
    def test_from_mp(self):
        """Test downloading structures from materialsproject.org."""
        # TODO Needs ensuring that the structure query gets the same
        # structure as we have downloaded.
        api_key = os.environ.get("MPI_KEY")

        for comp, species in self.TEST_SPECIES.items():
            with self.subTest(comp=comp):
                comp_file = os.path.join(files_dir, f"{comp}.txt")
                local_struct = SmactStructure.from_file(comp_file)
                mp_struct = SmactStructure.from_mp(species, api_key)
                self.assertEqual(local_struct, mp_struct)


class StructureDBTest(unittest.TestCase):
    """Test StructureDB interface."""

    TEST_DB = os.path.join(files_dir, "test_db.tmp")
    TEST_TABLE = "Structures"

    @classmethod
    def tearDownClass(cls):
        """Remove database files."""
        if os.path.exists(cls.TEST_DB):
            os.remove(cls.TEST_DB)

    def test_db_interface(self):
        """Test interfacing with database."""
        with self.subTest(msg="Instantiating database."):
            self.db = StructureDB(self.TEST_DB)

        with self.subTest(msg="Adding table."):
            self.assertTrue(self.db.add_table(self.TEST_TABLE))

        struct_file = os.path.join(files_dir, "CaTiO3.txt")
        struct = SmactStructure.from_file(struct_file)

        with self.subTest(msg="Adding structure to table."):
            self.assertTrue(self.db.add_struct(struct, self.TEST_TABLE))

        with self.subTest(msg="Getting structure from table."):
            struct_list = self.db.get_structs(struct.composition(), self.TEST_TABLE)
            self.assertEqual(len(struct_list), 1)
            self.assertEqual(struct_list[0], struct)

        struct_files = [os.path.join(files_dir, f"{x}.txt") for x in ["NaCl", "Fe"]]
        structs = [SmactStructure.from_file(fname) for fname in struct_files]

        with self.subTest(msg="Adding multiple structures to table."):
            self.assertTrue(self.db.add_structs(structs, self.TEST_TABLE))


class CationMutatorTest(unittest.TestCase):
    """Test the CationMutator class."""

    @classmethod
    def setUpClass(cls):
        """Set up the test initial structure and mutator."""
        cls.test_struct = SmactStructure.from_file(TEST_POSCAR)
        cls.test_mutator = CationMutator(lambda_json=TEST_LAMBDA_TAB)
        cls.test_pymatgen_mutator = CationMutator(lambda_json=None)

    def assertDataFrameEqual(self, df1: pd.DataFrame, df2: pd.DataFrame):
        """Assert that two pandas.DataFrames are equal.

        The indices and columns must also have identical labels.
        """
        self.assertTrue(df1.equals(df2))
        self.assertEqual(list(df1.index), list(df2.index))
        self.assertEqual(list(df1.columns), list(df2.columns))

    def test_lambda_tab_pop(self):
        """Test if lambda table is populated correctly."""
        lambda_dat = [
          [-5.0, 0.5, -5.0],
          [0.5, -5.0, 0.3],
          [-5.0, 0.3, -5.0], ]
        labels = ["A", "B", "C"]

        exp_lambda = pd.DataFrame(lambda_dat, index=labels, columns=labels)

        self.assertDataFrameEqual(self.test_mutator.lambda_tab, exp_lambda)

    def test_partition_func_Z(self):
        """Test the partition function for the whole table."""
        # 2e^0.5 + 2e^0.3 + 5e^{-5} \approx 6.0308499
        self.assertAlmostEqual(self.test_mutator.Z, 6.0308499)

    def test_pymatgen_lambda_import(self):
        """Test importing pymatgen lambda table."""
        self.assertIsInstance(self.test_pymatgen_mutator.lambda_tab, pd.DataFrame)

    def test_lambda_interface(self):
        """Test getting lambda values."""
        test_cases = [itertools.permutations(x) for x in [("A", "B"), ("A", "C"), ("B", "C")]]

        expected = [0.5, -5.0, 0.3]

        for test_case, expectation in zip(test_cases, expected):
            for spec_comb in test_case:
                s1, s2 = spec_comb
                with self.subTest(s1=s1, s2=s2):
                    self.assertEqual(self.test_mutator.get_lambda(s1, s2), expectation)

    def test_ion_mutation(self):
        """Test mutating an ion of a SmactStructure."""
        ca_file = os.path.join(files_dir, "CaTiO3.txt")
        ba_file = os.path.join(files_dir, "BaTiO3.txt")

        CaTiO3 = SmactStructure.from_file(ca_file)
        BaTiO3 = SmactStructure.from_file(ba_file)

        with self.subTest(s1="CaTiO3", s2="BaTiO3"):
            mutation = self.test_mutator._mutate_structure(CaTiO3, "Ca2+", "Ba2+")
            self.assertEqual(mutation, CaTiO3)

        na_file = os.path.join(files_dir, "NaCl.txt")
        NaCl = SmactStructure.from_file(na_file)

        with self.subTest(s1="Na1+Cl1-", s2="Na2+Cl1-"):
            with self.assertRaises(ValueError):
                self.test_mutator._mutate_structure(NaCl, "Na1+", "Na2+")

        # TODO Confirm functionality with more complex substitutions

    @unittest.skip("Not implemented.")
    def test_subs_probs(self):
        """Test acquiring multiple substitution probabilities."""
        ca_sub_probs = self.test_pymatgen_mutator.sub_probs("Ca2+")

    @unittest.skip("Not implemented.")
    def test_cond_probs(self):
        """Test acquiring multiple conditional substitution probabilities."""
        ca_sub_probs = self.test_pymatgen_mutator.cond_sub_probs("Ca2+")
        print(f"\n Ca sub prob: {ca_sub_probs['Ca2+']}")
        print(ca_sub_probs.describe())


if __name__ == "__main__":
    TestLoader = unittest.TestLoader()
    StructureTests = unittest.TestSuite()
    StructureTests.addTests(TestLoader.loadTestsFromTestCase(StructureTest))
    StructureTests.addTests(TestLoader.loadTestsFromTestCase(StructureDBTest))
    StructureTests.addTests(TestLoader.loadTestsFromTestCase(CationMutatorTest))
    runner = unittest.TextTestRunner()
    result = runner.run(StructureTests)