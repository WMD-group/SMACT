#!/usr/bin/env python

import unittest
import os
import pickle
import numpy as np

from smact import Species
from smact.builder import SmactStructure, StructureDB

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")


class StructureTest(unittest.TestCase):
    """`SmactStructure` testing."""

    TEST_STRUCT = os.path.join(files_dir, "test_struct")
    TEST_POSCAR = os.path.join(files_dir, "test_poscar.txt")

    TEST_SPECIES = {
        "CaTiO3": [('Ca', 2, 1), ('Ti', 4, 1), ('O', -2, 3)],
        "NaCl": [('Na', 1, 1), ('Cl', -1, 1)],
        "Fe": [('Fe', 0, 1)],
    }

    def test_as_poscar(self):
        """Test POSCAR generation."""
        for comp in self.TEST_SPECIES.keys():
            with self.subTest(comp=comp):
                comp_file = os.path.join(files_dir, f"{comp}.txt")
                with open(comp_file, "r") as f:
                    struct = SmactStructure.from_file(comp_file)
                    self.assertEqual(struct.as_poscar(), f.read())

    @staticmethod
    def gen_empty_struct(species):
        """Generate an empty set of arguments for `SmactStructure` testing."""
        lattice_mat = np.array([[0] * 3] * 3)

        if isinstance(species[0][0], str):
            species_strs = [
                "{ele}{charge}{sign}".format(
                    ele=spec[0],
                    charge=abs(spec[1]),
                    sign='+' if spec[1] >= 0 else '-',) for spec in species
            ]
        else:
            species_strs = [
                "{ele}{charge}{sign}".format(
                    ele=spec[0].symbol,
                    charge=abs(spec[0].oxidation),
                    sign='+' if spec[0].oxidation >= 0 else '-',) for spec in species
            ]

        sites = {spec: [[]] for spec in species_strs}
        return species, lattice_mat, sites

    def test_smactStruc_comp_key(self):
        """Test generation of a composition key for `SmactStructure`s."""
        s1 = SmactStructure(*self.gen_empty_struct([('Ba', 2, 2), ('O', -2, 1), ('F', -1, 2)]))
        s2 = SmactStructure(*self.gen_empty_struct([('Fe', 2, 1), ('Fe', 3, 2), ('O', -2, 4)]))

        Ba = Species('Ba', 2)
        O = Species('O', -2)
        F = Species('F', -1)
        Fe2 = Species('Fe', 2)
        Fe3 = Species('Fe', 3)

        s3 = SmactStructure(*self.gen_empty_struct([(Ba, 2), (O, 1), (F, 2)]))
        s4 = SmactStructure(*self.gen_empty_struct([(Fe2, 1), (Fe3, 2), (O, 4)]))

        Ba_2OF_2 = "Ba_2_2+F_2_1-O_1_2-"
        Fe_3O_4 = "Fe_2_3+Fe_1_2+O_4_2-"
        self.assertEqual(s1.composition(), Ba_2OF_2)
        self.assertEqual(s2.composition(), Fe_3O_4)
        self.assertEqual(s3.composition(), Ba_2OF_2)
        self.assertEqual(s4.composition(), Fe_3O_4)

    def test_smactStruc_from_file(self):
        """Test the `from_file` method of `SmactStructure`."""
        with open(self.TEST_STRUCT, 'rb') as f:
            s1 = pickle.load(f)

        s2 = SmactStructure.from_file(self.TEST_POSCAR)

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
        s1 = SmactStructure(*self.gen_empty_struct([('Fe', 2, 1), ('Fe', 3, 2), ('O', -2, 4)]))
        s1_stoics = {'Fe': 3, 'O': 4}
        s2 = SmactStructure(*self.gen_empty_struct([('Ba', 2, 2), ('O', -2, 1), ('F', -1, 2)]))
        s2_stoics = {'Ba': 2, 'O': 1, 'F': 2}

        for test, expected in [(s1, s1_stoics), (s2, s2_stoics)]:
            with self.subTest(species=test.species):
                self.assertEqual(test._get_ele_stoics(), expected)

    @unittest.skipUnless(os.environ.get("MPI_KEY"), "requires MPI key to be set.")
    def test_from_mp(self):
        """Test downloading structures from materialsproject.org."""
        # TODO Needs ensuring that the structure query gets the same
        # structure as we have downloaded.
        for comp, species in self.TEST_SPECIES.items():
            with self.subTest(comp=comp):
                comp_file = os.path.join(files_dir, f"{comp}.txt")
                with open(comp_file, "r") as f:
                    local_struct = SmactStructure.from_file(comp_file)
                    mp_struct = SmactStructure.from_mp(species)
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


if __name__ == "__main__":
    TestLoader = unittest.TestLoader()
    StructureTests = unittest.TestSuite()
    StructureTests.addTests(TestLoader.loadTestsFromTestCase(StructureTest))
    StructureTests.addTests(TestLoader.loadTestsFromTestCase(StructureDBTest))
    runner = unittest.TextTestRunner()
    result = runner.run(StructureTests)
