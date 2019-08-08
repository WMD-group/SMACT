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
import typing
from typing import List, Tuple, Union, Optional
from operator import itemgetter
# First example: using ase

from ase.spacegroup import crystal
from smact.lattice import Lattice, Site
from smact import Species
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.bond_valence import BVAnalyzer

MPI_KEY = os.environ.get("MPI_KEY")


class SmactStructure:
    """SMACT implementation inspired by pymatgen Structure class."""

    def __init__(
        self,
        species: List[Union[Tuple[str, int, int], Tuple[Species, int]]],
        structure: Optional[str] = None,):
        """Initialize class with constituent species."""
        self.species = self._sanitise_species(species)

        if structure is None:
            with MPRester(MPI_KEY) as m:
                eles = self._get_ele_stoics()

                formula = "".join(f"{ele}{stoic}" for ele, stoic in eles.items())
                structs = m.query(
                    criteria={"reduced_cell_formula": formula},
                    properties=["structure"],)

                if len(structs) == 0:
                    raise ValueError(
                        "Could not find composition in Materials Project Database, "
                        "please supply a structure."
                    )

                self.struct = structs[0]['structure']  # Default to first found structure
                bva = BVAnalyzer()
                self.struct = bva.get_oxi_state_decorated_structure(self.struct)

        else:
            # TODO Validate input structure
            self.struct = structure

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
        specs = SmactStructure._sanitise_species(species)

        with MPRester(MPI_KEY) as m:
            eles = {}
            for spec in specs:
                if spec[0] in eles:
                    eles[spec[0]] += spec[1]
                else:
                    eles[spec[0]] = spec[1]

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
            bva = BVAnalyzer()
            struct = bva.get_oxi_state_decorated_structure(struct)

    def _format_style(self, template: str, delim: Optional[str] = " "):
        """Format a given template string with the composition."""
        return delim.join(
            template.format(
                ele=specie[0],
                stoic=specie[1],
                charge=abs(specie[2]),
                sign='+' if specie[2] >= 0 else '-',) for specie in self.species
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
        return self._format_style(comp_style, delim="")

    def as_poscar(self):
        """Represent the structure as a POSCAR file compatible with VASP5."""
        poscar = self._format_style("{ele}{charge}{sign}") + "\n"

        poscar += "1.0\n"  # TODO Use actual lattice parameter

        poscar += "\n".join(" ".join(map(str, vec)) for vec in self.struct.lattice.matrix.tolist()) + "\n"

        poscar_specs = {}
        for spec in self.struct.species:
            spec_str = str(spec)
            if spec_str in poscar_specs:
                poscar_specs[spec_str] += 1
            else:
                poscar_specs[spec_str] = 1

        poscar += self._format_style("{ele}") + "\n"

        species_strs = self._format_style("{ele}{charge}{sign}")
        poscar += " ".join(str(poscar_specs[spec]) for spec in species_strs.split(" ")) + "\n"

        poscar += "Coordinates\n"
        for spec in species_strs.split(" "):
            for site in self.struct.sites:
                if site.species_string == spec:
                    poscar += " ".join(map(str, site.coords.tolist()))
                    poscar += f" {spec}\n"

        return poscar


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
