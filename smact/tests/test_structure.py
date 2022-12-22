"""Test structure prediction module"""

import itertools
import json
import logging
import os
import pickle
import unittest
from contextlib import contextmanager
from operator import itemgetter
from random import sample

import numpy as np
import pandas as pd
import pymatgen
from pandas.testing import assert_frame_equal, assert_series_equal
from pymatgen.analysis.structure_prediction.substitution_probability import (
    SubstitutionProbability,
)

import smact
from smact import Species
from smact.structure_prediction.database import StructureDB
from smact.structure_prediction.mutation import CationMutator
from smact.structure_prediction.prediction import StructurePredictor
from smact.structure_prediction.structure import SmactStructure

files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_STRUCT = os.path.join(files_dir, "test_struct")
TEST_POSCAR = os.path.join(files_dir, "test_poscar.txt")
TEST_PY_STRUCT = os.path.join(files_dir, "pymatgen_structure.json")
TEST_LAMBDA_JSON = os.path.join(files_dir, "test_lambda_tab.json")
TEST_LAMBDA_CSV = os.path.join(files_dir, "test_lambda_tab.csv")
TEST_PREDICTOR_DB = os.path.join(files_dir, "test_predictor.db")
TEST_PREDICTOR_TABLE = "TEST"


def generate_test_structure(comp: str) -> bool:
    """Generate a pickled test structure for comparison."""
    poscar_file = os.path.join(files_dir, f"{comp}.txt")
    s = SmactStructure.from_file(poscar_file)

    with open(TEST_STRUCT, "wb") as f:
        pickle.dump(s, f)

    with open(TEST_POSCAR, "w") as f:
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
        "CaTiO3": [("Ca", 2, 1), ("Ti", 4, 1), ("O", -2, 3)],
        "NaCl": [("Na", 1, 1), ("Cl", -1, 1)],
        "Fe": [("Fe", 0, 1)],
    }

    def assertStructAlmostEqual(
        self, s1: SmactStructure, s2: SmactStructure, places: int = 7
    ):
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
                with open(comp_file) as f:
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
                    sign="+" if spec[1] >= 0 else "-",
                )
                for spec in species
            ]
        else:
            species_strs = [
                "{ele}{charge}{sign}".format(
                    ele=spec[0].symbol,
                    charge=abs(spec[0].oxidation),
                    sign="+" if spec[0].oxidation >= 0 else "-",
                )
                for spec in species
            ]

        sites = {spec: [[]] for spec in species_strs}
        return species, lattice_mat, sites

    def test_from_py_struct(self):
        """Test generation of SmactStructure from a pymatgen Structure."""
        with open(TEST_PY_STRUCT) as f:
            d = json.load(f)
            py_structure = pymatgen.core.Structure.from_dict(d)

        with ignore_warnings(smact.structure_prediction.logger):
            s1 = SmactStructure.from_py_struct(py_structure)

        s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))

        self.assertStructAlmostEqual(s1, s2)

    def test_from_py_struct_icsd(self):
        """Test generation of SmactStructure from a pymatgen Structure using ICSD statistics to determine oxidation states."""
        with open(TEST_PY_STRUCT) as f:
            d = json.load(f)
            py_structure = pymatgen.core.Structure.from_dict(d)

        with ignore_warnings(smact.structure_prediction.logger):
            s1 = SmactStructure.from_py_struct(
                py_structure, determine_oxi="comp_ICSD"
            )

        s2 = SmactStructure.from_file(os.path.join(files_dir, "CaTiO3.txt"))

        self.assertStructAlmostEqual(s1, s2)

    def test_has_species(self):
        """Test determining whether a species is in a `SmactStructure`."""
        s1 = SmactStructure(
            *self._gen_empty_structure(
                [("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]
            )
        )

        self.assertTrue(s1.has_species(("Ba", 2)))
        self.assertFalse(s1.has_species(("Ba", 3)))
        self.assertFalse(s1.has_species(("Ca", 2)))

    def test_smactStruc_comp_key(self):
        """Test generation of a composition key for `SmactStructure`s."""
        s1 = SmactStructure(
            *self._gen_empty_structure(
                [("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]
            )
        )
        s2 = SmactStructure(
            *self._gen_empty_structure(
                [("Fe", 2, 1), ("Fe", 3, 2), ("O", -2, 4)]
            )
        )

        Ba = Species("Ba", 2)
        O = Species("O", -2)
        F = Species("F", -1)
        Fe2 = Species("Fe", 2)
        Fe3 = Species("Fe", 3)

        s3 = SmactStructure(
            *self._gen_empty_structure([(Ba, 2), (O, 1), (F, 2)])
        )
        s4 = SmactStructure(
            *self._gen_empty_structure([(Fe2, 1), (Fe3, 2), (O, 4)])
        )

        Ba_2OF_2 = "Ba_2_2+F_2_1-O_1_2-"
        Fe_3O_4 = "Fe_2_3+Fe_1_2+O_4_2-"
        self.assertEqual(s1.composition(), Ba_2OF_2)
        self.assertEqual(s2.composition(), Fe_3O_4)
        self.assertEqual(s3.composition(), Ba_2OF_2)
        self.assertEqual(s4.composition(), Fe_3O_4)

    def test_smactStruc_from_file(self):
        """Test the `from_file` method of `SmactStructure`."""
        with open(TEST_STRUCT, "rb") as f:
            s1 = pickle.load(f)

        s2 = SmactStructure.from_file(TEST_POSCAR)

        self.assertEqual(s1, s2)

    def test_equality(self):
        """Test equality determination of `SmactStructure`."""
        struct_files = [
            os.path.join(files_dir, f"{x}.txt") for x in ["CaTiO3", "NaCl"]
        ]
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
        s1 = SmactStructure(
            *self._gen_empty_structure(
                [("Fe", 2, 1), ("Fe", 3, 2), ("O", -2, 4)]
            )
        )
        s1_stoics = {"Fe": 3, "O": 4}
        s2 = SmactStructure(
            *self._gen_empty_structure(
                [("Ba", 2, 2), ("O", -2, 1), ("F", -1, 2)]
            )
        )
        s2_stoics = {"Ba": 2, "O": 1, "F": 2}

        for test, expected in [(s1, s1_stoics), (s2, s2_stoics)]:
            with self.subTest(species=test.species):
                self.assertEqual(
                    SmactStructure._get_ele_stoics(test.species), expected
                )

    @unittest.skipUnless(
        os.environ.get("MPI_KEY"), "requires MPI key to be set."
    )
    def test_from_mp(self):
        """Test downloading structures from materialsproject.org."""
        # TODO Needs ensuring that the structure query gets the same
        # structure as we have downloaded.
        # Need to modify the test for both legacy and next-gen queries
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
    TEST_MP_TABLE = "Structures1"

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
            try:
                self.db.add_table(self.TEST_TABLE)
            except Exception as e:
                self.fail(e)

        struct_file = os.path.join(files_dir, "CaTiO3.txt")
        struct = SmactStructure.from_file(struct_file)

        with self.subTest(msg="Adding structure to table."):
            try:
                self.db.add_struct(struct, self.TEST_TABLE)
            except Exception as e:
                self.fail(e)

        with self.subTest(msg="Getting structure from table."):
            struct_list = self.db.get_structs(
                struct.composition(), self.TEST_TABLE
            )
            self.assertEqual(len(struct_list), 1)
            self.assertEqual(struct_list[0], struct)

        struct_files = [
            os.path.join(files_dir, f"{x}.txt") for x in ["NaCl", "Fe"]
        ]
        structs = [SmactStructure.from_file(fname) for fname in struct_files]

        with self.subTest(msg="Adding multiple structures to table."):
            try:
                self.db.add_structs(structs, self.TEST_TABLE)
            except Exception as e:
                self.fail(e)

        test_with_species_args = [
            [("Na", 1)],
            [("Cl", -1)],
            [("Na", 1), ("Cl", -1)],
            [("Cl", -1), ("Na", 1)],
            [("Cl", -1)],
            [("Na", 1), ("Cl", 1)],
            [("O", -2)],
            [("Ca", 2), ("Ti", 4), ("O", -2)],
        ]

        test_with_species_exp = [
            [structs[0]],
            [structs[0]],
            [structs[0]],
            [structs[0]],
            [structs[0]],
            [],
            [struct],
            [struct],
        ]

        for spec, expected in zip(
            test_with_species_args, test_with_species_exp
        ):
            with self.subTest(msg=f"Retrieving species with {spec}"):
                self.assertEqual(
                    self.db.get_with_species(spec, self.TEST_TABLE), expected
                )

        mp_strucs = [
            pymatgen.core.Structure.from_file(os.path.join(files_dir, f))
            for f in ["CaTiO3.json", "NaCl.json", "Fe3O4.json"]
        ]
        mp_data = [
            {"material_id": mpid, "structure": s}
            for mpid, s in zip(["mp-4019", "mp-22862", "mp-19306"], mp_strucs)
        ]

        with self.subTest(msg="Testing adding downloaded MP structures."):
            added: int = self.db.add_mp_icsd(self.TEST_MP_TABLE, mp_data)
            self.assertEqual(added, 3)


class CationMutatorTest(unittest.TestCase):
    """Test the CationMutator class."""

    @classmethod
    def setUpClass(cls):
        """Set up the test initial structure and mutator."""
        cls.test_struct = SmactStructure.from_file(TEST_POSCAR)

        cls.test_mutator = CationMutator.from_json(
            lambda_json=TEST_LAMBDA_JSON
        )
        cls.test_pymatgen_mutator = CationMutator.from_json(
            lambda_json=None, alpha=lambda x, y: -5
        )

        # 5 random test species -> 5! test pairs
        cls.test_species = sample(cls.test_pymatgen_mutator.specs, 5)
        cls.test_pairs = list(
            itertools.combinations_with_replacement(cls.test_species, 2)
        )

        cls.pymatgen_sp = SubstitutionProbability(lambda_table=None, alpha=-5)

    def test_lambda_tab_pop(self):
        """Test if lambda table is populated correctly."""
        lambda_dat = [
            [-5.0, 0.5, -5.0],
            [0.5, -5.0, 0.3],
            [-5.0, 0.3, -5.0],
        ]
        labels = ["A", "B", "C"]

        exp_lambda = pd.DataFrame(lambda_dat, index=labels, columns=labels)

        assert_frame_equal(
            self.test_mutator.lambda_tab,
            exp_lambda,
            check_names=False,
        )

    def test_partition_func_Z(self):
        """Test the partition function for the whole table."""
        # 2e^0.5 + 2e^0.3 + 5e^{-5} \approx 6.0308499
        self.assertAlmostEqual(self.test_mutator.Z, 6.0308499)

    def test_pymatgen_lambda_import(self):
        """Test importing pymatgen lambda table."""
        self.assertIsInstance(
            self.test_pymatgen_mutator.lambda_tab, pd.DataFrame
        )

    def test_lambda_interface(self):
        """Test getting lambda values."""
        test_cases = [
            itertools.permutations(x)
            for x in [("A", "B"), ("A", "C"), ("B", "C")]
        ]

        expected = [0.5, -5.0, 0.3]

        for test_case, expectation in zip(test_cases, expected):
            for spec_comb in test_case:
                s1, s2 = spec_comb
                with self.subTest(s1=s1, s2=s2):
                    self.assertEqual(
                        self.test_mutator.get_lambda(s1, s2), expectation
                    )

    def test_ion_mutation(self):
        """Test mutating an ion of a SmactStructure."""
        ca_file = os.path.join(files_dir, "CaTiO3.txt")
        ba_file = os.path.join(files_dir, "BaTiO3.txt")

        CaTiO3 = SmactStructure.from_file(ca_file)
        BaTiO3 = SmactStructure.from_file(ba_file)

        with self.subTest(s1="CaTiO3", s2="BaTiO3"):
            mutation = self.test_mutator._mutate_structure(
                CaTiO3, "Ca2+", "Ba2+"
            )
            self.assertEqual(mutation, BaTiO3)

        na_file = os.path.join(files_dir, "NaCl.txt")
        NaCl = SmactStructure.from_file(na_file)

        with self.subTest(s1="Na1+Cl1-", s2="Na2+Cl1-"):
            with self.assertRaises(ValueError):
                self.test_mutator._mutate_structure(NaCl, "Na1+", "Na2+")

        # TODO Confirm functionality with more complex substitutions

    def test_sub_prob(self):
        """Test determining substitution probabilities."""
        for s1, s2 in self.test_pairs:
            with self.subTest(s1=s1, s2=s2):
                self.assertAlmostEqual(
                    self.pymatgen_sp.prob(s1, s2),
                    self.test_pymatgen_mutator.sub_prob(s1, s2),
                )

    def test_cond_sub_probs(self):
        """Test determining conditional substitution probabilities for a row."""
        for s1 in ["A", "B", "C"]:
            with self.subTest(s=s1):
                cond_sub_probs_test = self.test_mutator.cond_sub_probs(s1)

                vals = [
                    (s1, s2, self.test_mutator.cond_sub_prob(s1, s2))
                    for s2 in ["A", "B", "C"]
                ]

                test_df = pd.DataFrame(vals)
                test_df: pd.DataFrame = test_df.pivot(
                    index=0, columns=1, values=2
                )

                # Slice to convert to series
                assert_series_equal(cond_sub_probs_test, test_df.iloc[0])

    def test_cond_sub_prob(self):
        """Test determining conditional substitution probabilities."""
        for s1, s2 in self.test_pairs:
            with self.subTest(s1=s1, s2=s2):
                self.assertAlmostEqual(
                    self.pymatgen_sp.cond_prob(s1, s2),
                    self.test_pymatgen_mutator.cond_sub_prob(s1, s2),
                )

    def test_pair_corr(self):
        """Test determining conditional substitution probabilities."""
        for s1, s2 in self.test_pairs:
            with self.subTest(s1=s1, s2=s2):
                self.assertAlmostEqual(
                    self.pymatgen_sp.pair_corr(s1, s2),
                    self.test_pymatgen_mutator.pair_corr(s1, s2),
                )

    def test_from_df(self):
        """Test creating a CationMutator from an existing DataFrame."""
        lambda_df = pd.read_csv(TEST_LAMBDA_CSV, index_col=0)

        csv_test = CationMutator(lambda_df=lambda_df)
        assert_frame_equal(
            csv_test.lambda_tab,
            self.test_mutator.lambda_tab,
            check_names=False,
        )

    def test_complete_cond_probs(self):
        """Test getting all conditional probabilities."""
        pairs = itertools.product(["A", "B", "C"], repeat=2)

        vals = [
            (
                s1,
                s2,
                self.test_mutator.cond_sub_prob(s1, s2),
            )
            for s1, s2 in pairs
        ]

        cond_probs = pd.DataFrame(vals)
        cond_probs = cond_probs.pivot(index=0, columns=1, values=2)

        assert_frame_equal(self.test_mutator.complete_cond_probs(), cond_probs)

    def test_complete_sub_probs(self):
        """Test getting all probabilities."""
        pairs = itertools.product(["A", "B", "C"], repeat=2)

        vals = [
            (
                s1,
                s2,
                self.test_mutator.sub_prob(s1, s2),
            )
            for s1, s2 in pairs
        ]

        sub_probs = pd.DataFrame(vals)
        sub_probs = sub_probs.pivot(index=0, columns=1, values=2)

        assert_frame_equal(self.test_mutator.complete_sub_probs(), sub_probs)

    def test_complete_pair_corrs(self):
        """Test getting all pair correlations."""
        pairs = itertools.product(["A", "B", "C"], repeat=2)

        vals = [
            (
                s1,
                s2,
                self.test_mutator.pair_corr(s1, s2),
            )
            for s1, s2 in pairs
        ]

        pair_corrs = pd.DataFrame(vals)
        pair_corrs = pair_corrs.pivot(index=0, columns=1, values=2)

        assert_frame_equal(self.test_mutator.complete_pair_corrs(), pair_corrs)


class PredictorTest(unittest.TestCase):
    """Testing for the StructurePredictor wrapper."""

    @classmethod
    def setUpClass(cls):
        """Initialize the prerequisites."""
        cls.db = StructureDB(TEST_PREDICTOR_DB)
        # NOTE: This may break if the pymatgen lambda table is updated
        cls.cm = CationMutator.from_json()
        cls.table = TEST_PREDICTOR_TABLE

    def test_prediction(self):
        """Test prediction pipeline."""
        with self.subTest(msg="Testing predictor setup"):
            sp = StructurePredictor(self.cm, self.db, self.table)

        test_specs = [("Ca", 2), ("Ti", 4), ("O", -2)]

        # At time of creation, results with thresh=0 and include_same=True are as follows:
        #
        # Ca_1_2+O_4_2-Ti_2_4+: p=1.0, from Ca_1_2+O_4_2-Ti_2_4+
        # Ca_1_2+O_3_2-Ti_1_4+: p=0.0088406421433475, from Ca_1_2+O_3_2-Si_1_4+
        # Ca_1_2+O_3_2-Ti_1_4+: p=0.05502122697602804, from Ca_1_2+Mo_1_4+O_3_2-
        # Ca_1_2+O_3_2-Ti_1_4+: p=0.0445515892332118, from Ca_1_2+O_3_2-V_1_4+
        # Ca_2_2+O_4_2-Ti_1_4+: p=0.0274412135813993, from Fe_2_2+O_4_2-Ti_1_4+

        with self.subTest(msg="Acquiring predictions"):
            try:
                predictions = list(
                    sp.predict_structs(
                        test_specs, thresh=0.02, include_same=False
                    )
                )
            except Exception as e:
                self.fail(e)

        expected_comps = [
            "Ca_2_2+O_8_2-Sn_3_4+",
            "Ca_1_2+Mo_1_4+O_3_2-",
            "Ca_1_2+O_3_2-V_1_4+",
            "Ca_2_2+O_8_2-V_3_4+",
            "Fe_2_2+O_4_2-Ti_1_4+",
        ]
        expected_probs = [
            0.05604346851066936,
            0.05502122697602808,
            0.04455158923321182,
            0.04455158923321182,
            0.02744121358139931,
        ]

        # Sort from highest to lowest probability
        predictions.sort(key=itemgetter(1), reverse=True)

        parents = [x[2].composition() for x in predictions]

        with self.subTest(msg="Ensuring same parents"):
            self.assertEqual(parents, expected_comps)

        probs = [x[1] for x in predictions]

        with self.subTest(msg="Ensuring similar probabilities"):
            for p1, p2 in zip(probs, expected_probs):
                with self.subTest(p1=p1, p2=p2):
                    self.assertAlmostEqual(p1, p2)
