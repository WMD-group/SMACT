#!/usr/bin/env python
# Using the ase spacegroup module this can build the structure, from
# the composition, as defined in the smact_lattice module.
#TODO:
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

import os
import json
import re
import sqlite3
from contextlib import ContextDecorator
import typing
from typing import List, Tuple, Union, Optional, Dict, Sequence
from operator import itemgetter
# First example: using ase

from ase.spacegroup import crystal
from smact.lattice import Lattice, Site
from smact import Species
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.bond_valence import BVAnalyzer
import numpy as np


class SmactStructure:
    """SMACT implementation inspired by pymatgen Structure class."""

    def __init__(
        self,
        species: List[Union[Tuple[str, int, int], Tuple[Species, int]]],
        lattice_mat: np.ndarray,
        sites: Dict[str, List[List[float]]],
        lattice_param: Optional[float] = 1.0,
    ):
        """Initialize class with constituent species."""

        self.species = self._sanitise_species(species)
        self.lattice_mat = lattice_mat
        species_strs = self._format_style("{ele}{charge}{sign}")
        self.sites = {spec: sites[spec] for spec in species_strs.split(" ")}  # Sort sites
        self.lattice_param = lattice_param

    def __repr__(self):
        return self.as_poscar()

    @staticmethod
    def _sanitise_species(species: List[Union[Tuple[str, int, int], Tuple[Species, int]]]):
        """Sanitise and format a list of species."""

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
            sanit_species = [(x[0], x[2], x[1]) for x in species]  # Rearrange

        elif isinstance(species[0][0], Species):  # Species class variation of instantiation
            species.sort(key=lambda x: (x[0].symbol, -x[0].oxidation))
            sanit_species = [(x[0].symbol, x[1], x[0].oxidation) for x in species]

        else:
            raise TypeError(species_error)

        return sanit_species

    @staticmethod
    def from_mp(species: List[Union[Tuple[str, int, int], Tuple[Species, int]]]):
        """Create a SmactStructure using the first Materials Project entry for a composition."""

        MPI_KEY = os.environ.get("MPI_KEY")

        with MPRester(MPI_KEY) as m:
            eles = {}
            for spec in species:
                if spec[0] in eles:
                    eles[spec[0]] += spec[2]
                else:
                    eles[spec[0]] = spec[2]

            formula = "".join(f"{ele}{stoic}" for ele, stoic in eles.items())
            structs = m.query(
                criteria={"reduced_cell_formula": formula},
                properties=["structure"],)

            if len(structs) == 0:
                raise ValueError(
                    "Could not find composition in Materials Project Database, "
                    "please supply a structure."
                )

            struct = structs[0]['structure']  # Default to first found structure
            if not 0 in (spec[1] for spec in species):  # If everything's charged
                bva = BVAnalyzer()
                struct = bva.get_oxi_state_decorated_structure(struct)

        lattice_mat = struct.lattice.matrix

        lattice_param = 1.0  # TODO Use actual lattice parameter

        sites = {}
        for site in struct.sites:
            site_type = site.species_string
            # Add charge magnitude, for cases of unit charge
            if all([
                site_type[-2] not in map(str, range(10)),
                site_type[-1] in ("+", "-"),]):
                site_type = site_type[:-1] + '1' + site_type[-1]

            if site_type in sites:
                sites[site_type].append(site.coords.tolist())
            else:
                sites[site_type] = [site.coords.tolist()]

        return SmactStructure(species, lattice_mat, sites, lattice_param)

    @staticmethod
    def from_file(fname: str):
        """Create `SmactStructure` from a POSCAR file."""
        with open(fname, 'r') as f:
            return SmactStructure.from_poscar(f.read())

    @staticmethod
    def from_poscar(poscar: str):
        """Create `SmactStructure` from a POSCAR string."""

        lines = poscar.split("\n")

        # Find stoichiometry
        total_specs = [int(x) for x in lines[6].split(" ")]
        total_spec_sum = sum(total_specs)
        total_specs = [x / total_spec_sum for x in total_specs]
        total_spec_min = min(total_specs)
        total_specs = [round(x / total_spec_min) for x in total_specs]

        species = []
        for spec_str, stoic in zip(lines[0].split(" "), total_specs):
            charge_match = re.search(r"\d", spec_str)

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

        sites = {}
        for line in lines[8:]:
            if not line:  # EOF
                break

            split_line = line.split(" ")
            coords = [float(x) for x in split_line[:3]]
            spec = split_line[-1]

            if spec in sites:
                sites[spec].append(coords)
            else:
                sites[spec] = [coords]

        return SmactStructure(species, lattice, sites, lattice_param)

    @staticmethod
    def __get_sign(charge):
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
        include_ground: Optional[bool] = False,):
        """Format a given template string with the composition."""

        if include_ground:
            return delim.join(
                template.format(
                    ele=specie[0],
                    stoic=specie[1],
                    charge=abs(specie[2]),
                    sign="+" if specie[2] >= 0 else "-",) for specie in self.species
            )

        return delim.join(
            template.format(
                ele=specie[0],
                stoic=specie[1],
                charge=abs(specie[2]) if specie[2] != 0 else "",
                sign=self.__get_sign(specie[2]),
            ) for specie in self.species
        )

    def _get_ele_stoics(self):
        """Get the number of each element type in the compound, irrespective of oxidation state."""
        eles = {}
        for specie in self.species:
            if specie[0] in eles:
                eles[specie[0]] += specie[1]
            else:
                eles[specie[0]] = specie[1]
        return eles

    def composition(self):
        """Generate a key that describes the composition."""
        comp_style = "{ele}_{stoic}_{charge}{sign}"
        return self._format_style(comp_style, delim="", include_ground=True)

    def as_poscar(self):
        """Represent the structure as a POSCAR file compatible with VASP5."""
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
    """SQL Structure Database interface."""

    def __init__(self, db: str):
        """Set database name."""
        self.db = db

    def __enter__(self):
        """Initialize database connection."""
        self.conn = sqlite3.connect(self.db)
        self.cur = self.conn.cursor()

        return self.cur

    def __exit__(self, *args):
        """Close database connection."""
        self.conn.commit()
        self.conn.close()

    def add_table(self, table: str):
        """Add a table to the database."""
        with self as c:
            c.execute(f"""CREATE TABLE {table}
                (composition TEXT NOT NULL, structure TEXT NOT NULL)""",)
        return True

    def add_struct(self, struct: SmactStructure, table: str):
        """Add a `SmactStructure` to a table."""
        entry = (struct.composition(), struct.as_poscar())

        with self as c:
            c.execute(f"INSERT into {table} VALUES (?, ?)", entry)

        return True

    def add_structs(self, structs: Sequence[SmactStructure], table: str):
        """Add several `SmactStructure`s to a table."""
        with self as c:
            entries = [(
                struct.composition(),
                struct.as_poscar(),) for struct in structs]

            c.executemany(f"INSERT into {table} VALUES (?, ?)", entries)

        return True

    def get_structs(self, composition: str, table: str):
        """Get `SmactStructures` for a given composition."""
        with self as c:
            c.execute(
                f"SELECT structure FROM {table} WHERE composition = ?",
                (composition,),)
            structs = c.fetchall()
        return [SmactStructure.from_poscar(pos[0]) for pos in structs]


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
