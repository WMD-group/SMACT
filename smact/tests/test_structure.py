#!/usr/bin/env python

import unittest
import os
import numpy as np

from smact import Species
from smact.builder import SmactStructure, StructureDB


class StructureTest(unittest.TestCase):
    """Testing for the Structure predictor pipeline."""

    TEST_DB = "test.db"
    TEST_TABLE = "Structures"
    TEST_POSCAR = "test_poscar.test"

    @classmethod
    def tearDownClass(cls):
        """Remove temporary test files."""
        for fname in [cls.TEST_DB, cls.TEST_POSCAR]:
            if os.path.exists(fname):
                os.remove(fname)

    # ------------- SmactStructure --------------

    def test_as_poscar(self):
        """Test POSCAR generation."""
        test_cases = {
            "CaTiO3": [('Ca', 2, 1), ('Ti', 4, 1), ('O', -2, 3)],
            "NaCl": [('Na', 1, 1), ('Cl', -1, 1)],
            "Fe": [('Fe', 0, 1)],
        }
        for comp, species in test_cases.items():
            with self.subTest(comp=comp):
                with open(f"files/{comp}.txt", "r") as f:
                    struct = SmactStructure.from_mp(species)
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
        s1 = SmactStructure.from_mp([('Fe', 2, 1), ('Fe', 3, 2), ('O', -2, 4)])

        with open(self.TEST_POSCAR, 'w') as f:
            f.write(s1.as_poscar())

        s2 = SmactStructure.from_file(self.TEST_POSCAR)
        self.assertEqual(s1.species, s2.species)
        self.assertEqual(s1.lattice_mat.tolist(), s2.lattice_mat.tolist())
        self.assertEqual(s1.lattice_param, s2.lattice_param)
        self.assertDictEqual(s1.sites, s2.sites)

    # --------------- StructureDB ---------------

    def test_struct_addition(self):
        """Test adding to a database."""
        Db = StructureDB(self.TEST_DB)
        self.assertTrue(Db.add_table(self.TEST_TABLE))

        struct = SmactStructure.from_file("files/CaTiO3.txt")
        self.assertTrue(Db.add_struct(struct, self.TEST_TABLE))


if __name__ == "__main__":
    unittest.main()
