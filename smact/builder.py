#!/usr/bin/env python
# Using the ase spacegroup module this can build the structure, from
# the composition, as defined in the smact_lattice module.
# TODO:
# Estimate the cell parameters based on radii from tables.
# Add further types, Spinnel, Flourite, Delafossite ....
# Implement Structure class, c.f. dev_docs.
################################################################################
# Copyright Keith T Butler, Adam J Jackson    (2013)                           #
#                                                                              #
# This file is part of SMACT: builder.py is free software: you can             #
# redistribute it and/or modify it under the terms of the GNU General Public   #
# License as published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.                          #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for    #
# more details.                                                                #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
#                                                                              #
################################################################################

from collections import defaultdict
import itertools
import json
import math
import os
import re
import sqlite3
import logging
from operator import itemgetter
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from ase.spacegroup import crystal
import pandas as pd
import pymatgen
import pymatgen.analysis.structure_prediction as pymatgen_sp
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.ext.matproj import MPRester

from smact import Species
from smact.lattice import Lattice, Site

logger = logging.getLogger(__name__)


class SmactStructure:
    """SMACT implementation inspired by pymatgen Structure class.

    Handles basic structural and compositional information for a compound.
    Includes a lossless POSCAR-style specification for storing structures,
    allowing structures to be stored in files or databases, or to be pulled
    from the `Materials Project <https://www.materialsproject.org>`_.

    Attributes:
        species: A list of tuples describing the composition of the structure,
            stored as (element, oxidation, stoichiometry). The list is sorted
            alphabetically based on element symbol, and identical elements
            are sorted with highest charge first.
        lattice_mat: A numpy 3x3 array containing the lattice vectors.
        sites: A dictionary of {species: coords}, where species is a string
            representation of the species and coords is a list of position
            vectors, given as lists of length 3. For example:

            >>> s = SmactStructure.from_file('tests/files/NaCl.txt')
            >>> s.sites
            {'Cl1-': [[2.323624165, 1.643050405, 4.02463512]], 'Na1+': [[0.0, 0.0, 0.0]]}

        lattice_param: The lattice parameter.

    """

    def __init__(
      self,
      species: List[Union[Tuple[str, int, int], Tuple[Species, int]]],
      lattice_mat: np.ndarray,
      sites: Dict[str, List[List[float]]],
      lattice_param: Optional[float] = 1.0,
      sanitise_species: Optional[bool] = True,
    ):
        """Initialize structure with constituent species.

        Args:
            species: See :class:`~.SmactStructure`. May be supplied as either a list of
                (element, oxidation, stoichiometry) or (smact.Species, stoichiometry).
            lattice_mat: See :class:`~.SmactStructure`.
            sites: See :class:`~.SmactStructure`.
            lattice_param: See :class:`~.SmactStructure`.
            sanitise_species: Whether to sanitise check species. Should be `True` unless
                species have already been sanitised by a different constructor like
                :meth:`~.from_mp`.

        """
        self.species = self._sanitise_species(species) if sanitise_species else species

        self.lattice_mat = lattice_mat

        species_strs = self._format_style("{ele}{charge}{sign}")

        self.sites = {spec: sites[spec] for spec in species_strs.split(" ")}  # Sort sites

        self.lattice_param = lattice_param

    def __repr__(self):
        """Represent the structure as a POSCAR.

        Alias for :meth:`~.as_poscar`.
        """
        return self.as_poscar()

    def __eq__(self, other):
        """Determine equality of SmactStructures based on their attributes.

        :attr:`~.species`, :attr:`~.lattice_mat`, :attr:`~.lattice_param` and
        :attr:`~.sites` must all be equal for the comparison to be True.

        Note:
            For the SmactStructures to be equal their attributes must be
            *identical*. For example, it is insufficient that the two
            structures have the same space group or the same species;
            the site coordinates must be equal also.

        """
        if not isinstance(other, SmactStructure):
            return False
        return all([
          self.species == other.species,
          self.lattice_mat.tolist() == other.lattice_mat.tolist(),
          self.lattice_param == other.lattice_param,
          self.sites == other.sites,
          list(self.sites.keys()) == list(other.sites.keys()),
        ])

    @staticmethod
    def _sanitise_species(
      species: List[Union[Tuple[str, int, int], Tuple[Species, int]]],
    ) -> List[Tuple[str, int, int]]:
        """Sanitise and format a list of species.

        Args:
            species: See :meth:`~.__init__`.

        Returns:
            sanit_species: Sanity-checked species in the format of
            a list of (element, oxidation, stoichiometry).

        Raises:
            TypeError: species contains the wrong types.
            ValueError: species is either empty or contains tuples of
                incorrect length.

        """
        if not isinstance(species, list):
            raise TypeError(f"`species` must be a list, got {type(species)}.")
        if len(species) == 0:
            raise ValueError("`species` cannot be empty.")
        if not isinstance(species[0], tuple):
            raise TypeError(f"`species` must be a list of tuples, got list of {type(species[0])}.")

        species_error = (
          "`species` list of tuples must contain either "
          "2-tuples of Species objects and stoichiometries, "
          "or 3-tuples of elements, oxidations and stoichiometries."
        )
        if len(species[0]) not in {2, 3}:
            raise ValueError(species_error)

        if isinstance(species[0][0], str):  # String variation of instantiation
            species.sort(key=itemgetter(1), reverse=True)
            species.sort(key=itemgetter(0))
            sanit_species = species

        elif isinstance(species[0][0], Species):  # Species class variation of instantiation
            species.sort(key=lambda x: (x[0].symbol, -x[0].oxidation))
            sanit_species = [(x[0].symbol, x[0].oxidation, x[1]) for x in species]

        else:
            raise TypeError(species_error)

        return sanit_species

    @staticmethod
    def __parse_py_sites(
      structure: pymatgen.Structure,
    ) -> Tuple[Dict[str, List[List[float]]], List[Tuple[str, int, int]]]:
        """Parse the sites of a pymatgen Structure."""
        if not isinstance(structure, pymatgen.Structure):
            raise TypeError("structure must be a pymatgen.Structure instance.")

        sites = defaultdict(list)
        for site in structure.sites:
            site_type = site.species_string
            # Add charge magnitude, for cases of unit charge
            if all([
              site_type[-2] not in map(str, range(10)),
              site_type[-1] in ("+", "-"), ]):
                site_type = site_type[:-1] + '1' + site_type[-1]

            sites[site_type].append(site.coords.tolist())

        sites = dict(sites)

        # Find stoichiometry
        total_specs = [len(val) for val in sites.values()]
        total_spec_sum = sum(total_specs)
        total_specs = [x / total_spec_sum for x in total_specs]
        total_spec_min = min(total_specs)
        total_specs = [round(x / total_spec_min) for x in total_specs]

        species = []
        for spec, stoic in zip(sites.keys(), total_specs):
            charge_match = re.search(r"\d+", spec)

            if charge_match:
                charge_loc = charge_match.start()
                symb = spec[:charge_loc]
                charge = int(spec[-1] + spec[charge_loc:-1])
            else:
                symb = spec
                charge = 0

            species.append((symb, charge, stoic))

        return sites, species

    @staticmethod
    def from_py_struct(structure: pymatgen.Structure):
        """Create a SmactStructure from a pymatgen Structure object.

        Args:
            structure: A pymatgen Structure.

        Returns:
            :class:`~.SmactStructure`

        """
        if not isinstance(structure, pymatgen.Structure):
            raise TypeError("Structure must be a pymatgen.Structure instance.")

        try:
            bva = BVAnalyzer()
            struct = bva.get_oxi_state_decorated_structure(struct)
        except:
            logger.warn("Couldn't decorate structure with oxidation states.")

        sites, species = SmactStructure.__parse_py_sites(structure)

        lattice_mat = structure.lattice.matrix

        lattice_param = 1.0

        return SmactStructure(
          species,
          lattice_mat,
          sites,
          lattice_param,
          sanitise_species=True, )

    @staticmethod
    def from_mp(
      species: List[Union[Tuple[str, int, int], Tuple[Species, int]]],
      api_key: str, ):
        """Create a SmactStructure using the first Materials Project entry for a composition.

        Args:
            species: See :meth:`~.__init__`.
            api_key: A www.materialsproject.org API key.

        Returns:
            :class:`~.SmactStructure`

        """
        sanit_species = SmactStructure._sanitise_species(species)

        with MPRester(api_key) as m:
            eles = SmactStructure._get_ele_stoics(sanit_species)
            formula = "".join(f"{ele}{stoic}" for ele, stoic in eles.items())
            structs = m.query(
              criteria={"reduced_cell_formula": formula},
              properties=["structure"], )

            if len(structs) == 0:
                raise ValueError(
                  "Could not find composition in Materials Project Database, "
                  "please supply a structure."
                )

            struct = structs[0]['structure']  # Default to first found structure

        if 0 not in (spec[1] for spec in sanit_species):  # If everything's charged
            bva = BVAnalyzer()
            struct = bva.get_oxi_state_decorated_structure(struct)

        lattice_mat = struct.lattice.matrix

        lattice_param = 1.0  # TODO Use actual lattice parameter

        sites, _ = SmactStructure.__parse_py_sites(struct)

        return SmactStructure(
          sanit_species,
          lattice_mat,
          sites,
          lattice_param,
          sanitise_species=False, )

    @staticmethod
    def from_file(fname: str):
        """Create SmactStructure from a POSCAR file.

        Args:
            fname: The name of the POSCAR file.
                See :meth:`~.as_poscar` for format specification.

        Returns:
            :class:`~.SmactStructure`

        """
        with open(fname, 'r') as f:
            return SmactStructure.from_poscar(f.read())

    @staticmethod
    def from_poscar(poscar: str):
        """Create SmactStructure from a POSCAR string.

        Args:
            poscar: A SMACT-formatted POSCAR string.
                See :meth:`~.as_poscar` for format specification.

        Returns:
            :class:`~.SmactStructure`

        """
        lines = poscar.split("\n")

        # Find stoichiometry
        total_specs = [int(x) for x in lines[6].split(" ")]
        total_spec_sum = sum(total_specs)
        total_specs = [x / total_spec_sum for x in total_specs]
        total_spec_min = min(total_specs)
        total_specs = [round(x / total_spec_min) for x in total_specs]

        species = []
        for spec_str, stoic in zip(lines[0].split(" "), total_specs):
            charge_match = re.search(r"\d+", spec_str)

            if charge_match:
                charge_loc = charge_match.start()
                symb = spec_str[:charge_loc]
                charge = int(spec_str[-1] + spec_str[charge_loc:-1])
            else:
                symb = spec_str
                charge = 0

            species.append((symb, charge, stoic))

        lattice_param = float(lines[1])

        lattice = np.array([[float(point) for point in line.split(" ")] for line in lines[2:5]])

        sites = defaultdict(list)
        for line in lines[8:]:
            if not line:  # EOF
                break

            split_line = line.split(" ")
            coords = [float(x) for x in split_line[:3]]
            spec = split_line[-1]

            sites[spec].append(coords)

        sites = dict(sites)

        return SmactStructure(species, lattice, sites, lattice_param)

    @staticmethod
    def __get_sign(charge: int) -> str:
        """Get string representation of a number's sign.

        Args:
            charge: The number whose sign to derive.

        Returns:
            Sign; either '+', '-' or '' for neutral.

        """
        if charge > 0:
            return '+'
        elif charge < 0:
            return '-'
        else:
            return ''

    def _format_style(
      self,
      template: str,
      delim: Optional[str] = " ",
      include_ground: Optional[bool] = False, ) -> str:
        """Format a given template string with the composition.

        Formats a python template string with species information,
        with each species separated by a given delimiter.

        Args:
            template: Template string to format, using python's
                curly brackets notation. Supported keywords are
                `ele` for the elemental symbol, `stoic` for the
                stoichiometry, `charge` for the absolute value
                of oxidation state and `sign` for the 
                oxidation state's sign.
            delim: The delimeter between species' templates.
            include_ground: Whether to include the charge and sign
                of neutral species. See also :meth:`~.__get_sign`.
        
        Returns:
            String of templates formatted for each species, separated
                by `delim`.
        
        Examples:
            >>> s = SmactStructure.from_file('tests/files/CaTiO3.txt')
            >>> template = '{stoic}x{ele}{charge}{sign}'
            >>> print(s._format_style(template))
            1xCa2+ 3xO2- 1xTi4+

        """
        if include_ground:
            return delim.join(
              template.format(
                ele=specie[0],
                stoic=specie[2],
                charge=abs(specie[1]),
                sign="+" if specie[1] >= 0 else "-",
              ) for specie in self.species
            )

        return delim.join(
          template.format(
            ele=specie[0],
            stoic=specie[2],
            charge=abs(specie[1]) if specie[1] != 0 else "",
            sign=self.__get_sign(specie[1]),
          ) for specie in self.species
        )

    @staticmethod
    def _get_ele_stoics(species: List[Tuple[str, int, int]]) -> Dict[str, int]:
        """Get the number of each element type in the compound, irrespective of oxidation state.

        Args:
            species: See :meth:`~.__init__`.

        Returns:
            eles: Dictionary of {element: stoichiometry}.

        Examples:
            >>> species = [('Fe', 2, 1), ('Fe', 3, 2), ('O', -2, 4)]
            >>> print(SmactStructure._get_ele_stoics(species))
            {'Fe': 3, 'O': 4}

        """
        eles = defaultdict(int)
        for specie in species:
            eles[specie[0]] += specie[2]

        return dict(eles)

    def composition(self) -> str:
        """Generate a key that describes the composition.

        Key format is '{element}_{stoichiometry}_{charge}{sign}'
        with no delimiter, *sans brackets*. Species are ordered as stored within
        the structure, see :class:`~.SmactStructure`.

        Returns:
            Key describing constituent species.

        Examples:
            >>> s = SmactStructure.from_file('tests/files/CaTiO3.txt')
            >>> print(s.composition())
            Ca_1_2+O_3_2-Ti_1_4+

        """
        comp_style = "{ele}_{stoic}_{charge}{sign}"
        return self._format_style(comp_style, delim="", include_ground=True)

    def as_poscar(self) -> str:
        """Represent the structure as a POSCAR file compatible with VASP5.

        The POSCAR format adopted is as follows:

        The first line contains the species' names separated by a whitespace.
        The second through fourth line, inclusive, contain the lattice
        matrix: each line contains a lattice vector, with elements
        separated by a whitespace.
        The fifth line contains the elements' names separated by a whitespace.
        If more than one oxidation state exists for an element, the element
        appears multiple times; once for each oxidation state.
        The sixth line is the string 'Cartesian'.
        The seventh line onwards are the Cartesian coordinates of each site,
        separated by a whitespace. In addition, at the end of each line is the
        species' name, separated by a whitespace.

        For examples of this format, see the text files under tests/files.

        Returns:
            str: POSCAR-style representation of the structure.

        """
        poscar = self._format_style("{ele}{charge}{sign}") + "\n"

        poscar += f"{self.lattice_param}\n"

        poscar += "\n".join(" ".join(map(str, vec)) for vec in self.lattice_mat.tolist()) + "\n"

        spec_count = {spec: len(coords) for spec, coords in self.sites.items()}

        poscar += self._format_style("{ele}") + "\n"

        species_strs = self._format_style("{ele}{charge}{sign}")
        poscar += " ".join(str(spec_count[spec]) for spec in species_strs.split(" ")) + "\n"

        poscar += "Cartesian\n"
        for spec, coords in self.sites.items():
            for coord in coords:
                poscar += " ".join(map(str, coord))
                poscar += f" {spec}\n"

        return poscar


class StructureDB:
    """SQLite Structure Database interface.

    Acts as a context manager for database interfacing
    and wraps several useful SQLite commands within
    methods.

    Attributes:
        db: The database name.
        conn: The database connection. Only open when
            used as a context manager.
        cur: The database connection cursor. Only usable
            when class implemented as context manager.

    Examples:
        Connecting to a database in memory:
        >>> DB = StructureDB(':memory:')
        >>> with DB as c:
        ...     _ = c.execute("CREATE TABLE test (id, val)")
        ...     c.execute("SELECT * FROM test").fetchall()
        []
        >>> DB.cur.execute("SELECT * FROM test").fetchall()
        Traceback (most recent call last):
            ...
        sqlite3.ProgrammingError: Cannot operate on a closed database.

    """

    def __init__(self, db: str):
        """Set database name."""
        self.db = db

    def __enter__(self) -> sqlite3.Cursor:
        """Initialize database connection."""
        self.conn = sqlite3.connect(self.db)
        self.cur = self.conn.cursor()

        return self.cur

    def __exit__(self, *args):
        """Close database connection."""
        self.conn.commit()
        self.conn.close()

    def add_table(self, table: str) -> bool:
        """Add a table to the database.

        Args:
            table: The name of the table to add

        Returns:
            bool: Whether the operation was successful.

        """
        with self as c:
            c.execute(
              f"""CREATE TABLE {table}
                (composition TEXT NOT NULL, structure TEXT NOT NULL)""",
            )
        return True

    def add_struct(self, struct: SmactStructure, table: str) -> bool:
        """Add a SmactStructure to a table.

        Args:
            struct: The :class:`~.SmactStructure` to add.
            table: The name of the table to add the structure to.

        Returns:
            bool: Whether the operation was successful.

        """
        entry = (struct.composition(), struct.as_poscar())

        with self as c:
            c.execute(f"INSERT into {table} VALUES (?, ?)", entry)

        return True

    def add_structs(self, structs: Sequence[SmactStructure], table: str) -> bool:
        """Add several SmactStructures to a table.

        Args:
            structs: Iterable of :class:`~.SmactStructure`s to add to table.
            table: The name of the table to add the structs to.

        Returns:
            bool: Whether the operation was successful.

        """
        with self as c:
            entries = [(
              struct.composition(),
              struct.as_poscar(), ) for struct in structs]

            c.executemany(f"INSERT into {table} VALUES (?, ?)", entries)

        return True

    def get_structs(self, composition: str, table: str) -> List[SmactStructure]:
        """Get SmactStructures for a given composition.

        Args:
            composition: The composition to search for.
                See :meth:`SmactStructure.composition`.
            table: The name of the table in which to search.

        Returns:
            A list of :class:`~.SmactStructure` s.

        """
        with self as c:
            c.execute(
              f"SELECT structure FROM {table} WHERE composition = ?",
              (composition, ), )
            structs = c.fetchall()
        return [SmactStructure.from_poscar(pos[0]) for pos in structs]


class CationMutator:
    """Handles cation mutation of SmactStructures based on substitution probability.
    
    Based on the algorithm presented in::
        Hautier, G., Fischer, C., Ehrlacher, V., Jain, A., and Ceder, G. (2011)
        Data Mined Ionic Substitutions for the Discovery of New Compounds.
        Inorganic Chemistry, 50(2), 656-663.
        `doi:10.1021/ic102031h <https://pubs.acs.org/doi/10.1021/ic102031h>`_
    
    """

    def __init__(
      self,
      init_structure: SmactStructure,
      lambda_json: Optional[str] = None,
      alpha: Optional[Callable[[str, str], float]] = (lambda s1, s2: -5.0),
    ):
        """Assign attributes and get lambda table."""
        if lambda_json is not None:
            with open(lambda_json, 'r') as f:
                lambda_dat = json.load(f)
        else:
            # Get pymatgen lambda table
            py_sp_dir = os.path.dirname(pymatgen_sp.__file__)
            pymatgen_lambda = os.path.join(py_sp_dir, "data", "lambda.json")
            with open(pymatgen_lambda, 'r') as f:
                lambda_dat = json.load(f)

        # Convert lambda table to pandas DataFrame
        lambda_dat = [tuple(x) for x in lambda_dat]
        lambda_df = pd.DataFrame(lambda_dat)

        self.lambda_tab = lambda_df.pivot(index=0, columns=1, values=2)

        self.init_structure = init_structure
        self.alpha = alpha

        # Sum over all unique pair values
        specs = set(
          itertools.chain.from_iterable(
            set(getattr(self.lambda_tab, x)) for x in ["columns", "index"]
          )
        )
        pairs = itertools.combinations_with_replacement(specs, 2)
        self.Z = sum(math.exp(self.get_lambda(*pair)) for pair in pairs)

    def _species_in_lambda(self, species: str) -> bool:
        """Determine whether a species has records in the lambda table."""
        return any(species in getattr(self.lambda_tab, x) for x in ["columns", "index"])

    def get_lambda(self, s1: str, s2: str) -> float:
        """Get lambda values corresponding to a pair of species."""
        l = np.nan

        if all(self._species_in_lambda(s) for s in [s1, s2]):
            # s1 and s2 exist in lambda table
            for si, sj in itertools.permutations([s1, s2]):
                try:
                    l = self.lambda_tab.at[si, sj]
                except KeyError:
                    continue

                if not np.isnan(l):
                    break

        if np.isnan(l):
            l = self.alpha(s1, s2)

        return l

    @staticmethod
    def _mutate_structure(
      structure: SmactStructure,
      init_species: str,
      final_species: str, ) -> SmactStructure:
        """Mutate an ion within a SmactStructure."""

        def parse_spec(species: str) -> Tuple[str, int]:
            """Parse a species string into its element and charge."""
            ele = re.match(r"[A-Za-z]+", species).group(0)

            charge_match = re.search(r"\d+", species)
            charge = int(charge_match.group(0)) if charge_match else 0

            if '-' in species:
                charge *= -1

            return ele, charge

        init_spec_tup = parse_spec(init_species)
        struct_spec_tups = [spec[:2] for spec in structure.species]
        spec_loc = struct_spec_tups.index(init_spec_tup)

        final_spec_tup = parse_spec(final_species)

        # # Resolve stoichiometries to maintain charge neutrality

        ## Nonsensical without a change in sites / structure!

        # charge_ratio = final_spec_tup[1] / init_spec_tup[1]
        # stoics = [spec[2] for spec in structure.species]
        # stoics = [x * charge_ratio for x in stoics]

        # stoic_sum = sum(stoics)
        # stoics = [x / stoic_sum for x in stoics]
        # min_stoic = min(stoics)
        # stoics = [round(x / min_stoic) for x in stoics]

        # structure.species = [(*x[:2], y) for x, y in zip(structure.species, stoics)]

        # Replace species tuple
        structure.species[spec_loc] = (*final_spec_tup, structure.species[spec_loc][2])

        # Check for charge neutrality
        if sum(x[1] * x[2] for x in structure.species) != 0:
            raise ValueError("New structure is not charge neutral.")

        # Sort species again
        structure.species.sort(key=itemgetter(1), reverse=True)
        structure.species.sort(key=itemgetter(0))

        # Replace sites
        structure.sites[final_species] = structure.sites.pop(init_species)
        # And sort
        species_strs = structure._format_style("{ele}{charge}{sign}").split(" ")
        structure.sites = {spec: structure.sites[spec] for spec in species_strs}

        return structure

    def unary_substitute(self, structure: SmactStructure) -> Tuple[SmactStructure, float]:
        """Substitute a single ion within a structure."""
        pass


def cubic_perovskite(species, cell_par=[6, 6, 6, 90, 90, 90], repetitions=[1, 1, 1]):
    """
    Build a perovskite cell using the crystal function in ASE.

    Args:
        species (str): Element symbols
        cell_par (list): Six floats/ints specifying 3 unit cell lengths and 3 unit cell angles.
        repetitions (list): Three floats specifying the expansion of the cell in x,y,z directions.
    Returns:
        SMACT Lattice object of the unit cell,
        ASE crystal system of the unit cell.

    """
    system = crystal((species),
                     basis=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0)],
                     spacegroup=221,
                     size=repetitions,
                     cellpar=cell_par)

    sites_list = []
    oxidation_states = [[2]] + [[4]] + [[-2]] * 3
    for site in zip(system.get_scaled_positions(), oxidation_states):
        sites_list.append(Site(site[0], site[1]))

    return Lattice(sites_list, oxidation_states), system


def wurtzite(species, cell_par=[2, 2, 6, 90, 90, 120], repetitions=[1, 1, 1]):
    """
    Build a wurzite cell using the crystal function in ASE.

    Args:
        species (str): Element symbols
        cell_par (list): Six floats/ints specifying 3 unit cell lengths and 3 unit cell angles.
        repetitions (list): Three floats specifying the expansion of the cell in x,y,z directions.
    Returns:
        SMACT Lattice object of the unit cell,
        ASE crystal system of the unit cell.

    """
    system = crystal((species),
                     basis=[(2. / 3., 1. / 3., 0), (2. / 3., 1. / 3., 5. / 8.)],
                     spacegroup=186,
                     size=repetitions,
                     cellpar=[3, 3, 6, 90, 90, 120])

    sites_list = []
    oxidation_states = [[1], [2], [3], [4]] + [[-1], [-2], [-3], [-4]]

    for site in zip(system.get_scaled_positions(), oxidation_states):
        sites_list.append(Site(site[0], site[1]))
    return Lattice(sites_list, oxidation_states), system
