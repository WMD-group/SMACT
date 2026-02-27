"""Minimalist structure representation for comprehensible manipulation."""

from __future__ import annotations

import logging
import os
import re
from collections import defaultdict
from functools import reduce
from math import gcd
from operator import itemgetter
from typing import cast as _cast

import numpy as np
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core import SETTINGS
from pymatgen.core import Structure as pmg_Structure

try:
    from pymatgen.ext.matproj import MPRester

    HAS_LEGACY_MPRESTER = True
except ImportError:  # pragma: no cover
    HAS_LEGACY_MPRESTER = False
    MPRester = None

try:
    from mp_api.client import MPRester as MPResterNew

    HAS_MP_API = True
except ImportError:  # pragma: no cover
    HAS_MP_API = False
    MPResterNew = None
from pymatgen.transformations.standard_transformations import (
    OxidationStateDecorationTransformation,
)

import smact

from .utilities import get_sign

logger = logging.getLogger(__name__)


class SmactStructure:
    """
    SMACT implementation inspired by pymatgen Structure class.

    Handles basic structural and compositional information for a compound.
    Includes a lossless POSCAR-style specification for storing structures,
    allowing structures to be stored in files or databases, or to be pulled
    from the `Materials Project <https://www.materialsproject.org>`_.

    Attributes:
    ----------
        species: A list of tuples describing the composition of the structure,
            stored as (element, oxidation, stoichiometry). The list is sorted
            alphabetically based on element symbol, and identical elements
            are sorted with highest charge first.
        lattice_mat: A numpy 3x3 array containing the lattice vectors.
        sites: A dictionary of {species: coords}, where species is a string
            representation of the species and coords is a list of position
            vectors, given as lists of length 3. For example:

            >>> s = SmactStructure.from_file("tests/files/NaCl.txt")
            >>> s.sites
            {'Cl1-': [[2.323624165, 1.643050405, 4.02463512]], 'Na1+': [[0.0, 0.0, 0.0]]}

        lattice_param: The lattice parameter.

    """

    def __init__(
        self,
        species: list[tuple[str, int, int] | tuple[smact.Species, int]],
        lattice_mat: np.ndarray,
        sites: dict[str, list[list[float]]],
        lattice_param: float | None = 1.0,
        sanitise_species: bool | None = True,
    ) -> None:
        """
        Initialize structure with constituent species.

        Args:
        ----
            species: See :class:`~.SmactStructure`. May be supplied as either a list of
                (element, oxidation, stoichiometry) or (:class:`~smact.Species`, stoichiometry).
            lattice_mat: See :class:`~.SmactStructure`.
            sites: See :class:`~.SmactStructure`.
            lattice_param: See :class:`~.SmactStructure`.
            sanitise_species: Whether to sanitise check species. Should be `True` unless
                species have already been sanitised by a different constructor like
                :meth:`~.from_mp`.

        """
        # _sanitise_species always returns list[tuple[str, int, int]]; cast for type safety
        self.species: list[tuple[str, int, int]] = (
            self._sanitise_species(species) if sanitise_species else _cast("list[tuple[str, int, int]]", species)
        )

        self.lattice_mat = lattice_mat

        self.sites = {spec: sites[spec] for spec in self.get_spec_strs()}  # Sort sites

        self.lattice_param = lattice_param

    def __repr__(self) -> str:
        """
        Represent the structure as a POSCAR.

        Alias for :meth:`~.as_poscar`.

        """
        return self.as_poscar()

    def __eq__(self, other: object) -> bool:
        """
        Determine equality of SmactStructures based on their attributes.

        :attr:`~.species`, :attr:`~.lattice_mat`, :attr:`~.lattice_param` and
        :attr:`~.sites` must all be equal for the comparison to be True.

        Note:
        ----
            For the SmactStructures to be equal their attributes must be
            *identical*. For example, it is insufficient that the two
            structures have the same space group or the same species;
            the site coordinates must be equal also.

        """
        if not isinstance(other, SmactStructure):
            return False

        if list(self.sites.keys()) != list(other.sites.keys()):
            return False

        sites_equal = True
        for si, sj in zip(self.sites.values(), other.sites.values(), strict=True):
            if len(si) != len(sj):
                sites_equal = False
                break
            for ci, cj in zip(si, sj, strict=True):
                if not np.allclose(ci, cj, atol=1e-7, rtol=0):
                    sites_equal = False
                    break

        return all(
            [
                self.species == other.species,
                np.array_equal(self.lattice_mat, other.lattice_mat),
                self.lattice_param == other.lattice_param,
                sites_equal,
            ]
        )

    def __hash__(self) -> int:
        """
        Provide a hash consistent with :meth:`__eq__`.

        The hash is computed from immutable representations of the
        structure's key attributes. Coordinates are rounded to 1e-7 to
        match the tolerance used in equality checks (np.allclose with
        atol=1e-7). Note: mutating a hashed instance (changing sites,
        species, or lattice) will change its logical identity and can
        break dictionary/set invariants; prefer using immutable objects
        as dict keys.
        """
        # species is a list of tuples -> make it a tuple of tuples
        species_t = tuple(tuple(s) for s in self.species)

        # lattice_mat is a numpy array -> convert to tuple of tuples
        lattice_t = tuple(tuple(row) for row in self.lattice_mat.tolist())

        # sites: preserve insertion/order of keys (constructor enforces order)
        # Round coordinates to 1e-7 to match __eq__ tolerance
        sites_t = tuple(
            (
                spec,
                tuple(tuple(round(c, 7) for c in coord) for coord in coords),
            )
            for spec, coords in self.sites.items()
        )

        return hash((species_t, lattice_t, self.lattice_param, sites_t))

    @staticmethod
    def _sanitise_species(
        species: list[tuple[str, int, int] | tuple[smact.Species, int]],
    ) -> list[tuple[str, int, int]]:
        """
        Sanitise and format a list of species.

        Args:
        ----
            species: See :meth:`~.__init__`.

        Returns:
        -------
            sanit_species: Sanity-checked species in the format of
            a list of (element, oxidation, stoichiometry).

        Raises:
        ------
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
            str_species = _cast("list[tuple[str, int, int]]", species)
            str_species.sort(key=itemgetter(1), reverse=True)
            str_species.sort(key=itemgetter(0))
            sanit_species: list[tuple[str, int, int]] = str_species

        elif isinstance(species[0][0], smact.Species):  # Species class variation of instantiation
            sp_species = _cast("list[tuple[smact.Species, int]]", species)
            sp_species.sort(key=lambda x: (x[0].symbol, -x[0].oxidation))
            sanit_species = [(x[0].symbol, x[0].oxidation, x[1]) for x in sp_species]

        else:
            raise TypeError(species_error)

        return sanit_species

    @staticmethod
    def __parse_py_sites(
        structure: pmg_Structure,
    ) -> tuple[dict[str, list[list[float]]], list[tuple[str, int, int]]]:
        """
        Parse the sites of a pymatgen Structure.

        Args:
        ----
            structure: A :class:`pmg_Structure` instance.

        Returns:
        -------
            sites (dict): In which a key is a species string
                and its corresponding value is a list of the coordinates
                that species occupies in the supercell. The coordinates
                are represented by lists containing three elements: one
                for each spatial dimension.
            species (list): A list of each species in the structure,
                represented by a tuple of (element, charge, stoichiometry).

        """
        if not isinstance(structure, pmg_Structure):
            raise TypeError(f"Expected pymatgen.core.Structure, got {type(structure).__name__}")

        sites = defaultdict(list)
        for site in structure.sites:
            site_type = site.species_string
            # Add charge magnitude, for cases of unit charge
            if all(
                [
                    site_type[-2] not in map(str, range(10)),
                    site_type[-1] in ("+", "-"),
                ]
            ):
                site_type = site_type[:-1] + "1" + site_type[-1]

            sites[site_type].append(site.coords.tolist())

        sites = dict(sites)

        # Find stoichiometry
        total_specs = [len(val) for val in sites.values()]
        hcf = reduce(gcd, total_specs)
        total_specs = [int(x / hcf) for x in total_specs]

        species = []
        for spec, stoic in zip(sites.keys(), total_specs, strict=True):
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
    def from_py_struct(structure: pmg_Structure, determine_oxi: str = "BV") -> SmactStructure:
        """
        Create a SmactStructure from a pymatgen Structure object.

        Args:
        ----
            structure: A pymatgen Structure.
            determine_oxi (str): The method to determine the assignments oxidation states in the structure.
                Options are 'BV', 'comp_ICSD','both' for determining the oxidation states by bond valence,
                ICSD statistics or trial both sequentially, respectively.

        Returns:
        -------
            :class:`~.SmactStructure`

        """
        if not isinstance(structure, pmg_Structure):
            raise TypeError(f"Expected pymatgen.core.Structure, got {type(structure).__name__}")

        if determine_oxi == "BV":
            bva = BVAnalyzer()
            struct = bva.get_oxi_state_decorated_structure(structure)

        elif determine_oxi == "comp_ICSD":
            comp = structure.composition
            oxi_transform = OxidationStateDecorationTransformation(comp.oxi_state_guesses()[0])
            struct = oxi_transform.apply_transformation(structure)
            logger.info("Charge assigned based on ICSD statistics")

        elif determine_oxi == "both":
            try:
                bva = BVAnalyzer()
                struct = bva.get_oxi_state_decorated_structure(structure)
                logger.info("Oxidation states assigned using bond valence")
            except (ValueError, RuntimeError):
                comp = structure.composition
                oxi_transform = OxidationStateDecorationTransformation(comp.oxi_state_guesses()[0])
                struct = oxi_transform.apply_transformation(structure)
                logger.info("Oxidation states assigned based on ICSD statistics")
        elif determine_oxi == "predecorated":
            struct = structure

        else:
            raise ValueError(
                f"Argument for 'determine_oxi', <{determine_oxi}> is not valid. Choose either 'BV','comp_ICSD','both' or 'predecorated'."
            )

        sites, species = SmactStructure.__parse_py_sites(struct)

        lattice_mat = struct.lattice.matrix

        lattice_param = 1.0

        return SmactStructure(
            _cast("list[tuple[str, int, int] | tuple[smact.Species, int]]", species),
            lattice_mat,
            sites,
            lattice_param,
            sanitise_species=True,
        )

    @staticmethod
    def from_mp(
        species: list[tuple[str, int, int] | tuple[smact.Species, int]],
        api_key: str | None = None,
        determine_oxi: str = "BV",
    ) -> SmactStructure:
        """
        Create a SmactStructure using the first Materials Project entry for a composition.

        Args:
        ----
            species: See :meth:`~.__init__`.
            determine_oxi (str): The method to determine the assignments oxidation states in the structure. Options are 'BV', 'comp_ICSD','both' for determining the oxidation states by bond valence, ICSD statistics or trial both sequentially, respectively.
            api_key (str| None): A www.materialsproject.org API key.

        Returns:
        -------
            :class:`~.SmactStructure`

        """
        sanit_species = SmactStructure._sanitise_species(species)
        eles = SmactStructure._get_ele_stoics(sanit_species)
        formula = "".join(f"{ele}{stoic}" for ele, stoic in eles.items())

        if api_key is None:
            api_key = os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY")
        if api_key is None:
            raise ValueError(
                "No Materials Project API key found. Set the MP_API_KEY or PMG_MAPI_KEY environment variable."
            )

        # Legacy API routine
        if len(api_key) != 32:
            if not HAS_LEGACY_MPRESTER:
                raise ImportError(
                    "pymatgen legacy MPRester is not available. "
                    "Install pymatgen >= 2022.1 with legacy MP API support or use mp-api."
                )
            if MPRester is None:  # pragma: no cover
                raise ImportError("pymatgen legacy MPRester is not available.")
            with MPRester(api_key) as m:
                structs = m.query(
                    criteria={"reduced_cell_formula": formula},
                    properties=["structure"],
                )

        elif HAS_MP_API:
            # New API routine
            if MPResterNew is None:  # pragma: no cover
                raise ImportError("mp-api MPRester is not available.")
            with MPResterNew(api_key, use_document_model=False) as m:
                structs = m.materials.summary.search(formula=formula, fields=["structure"])
        else:
            if not HAS_LEGACY_MPRESTER:
                raise ImportError(
                    "Neither mp-api nor pymatgen legacy MPRester is available. Install mp-api: `pip install mp-api`."
                )
            if MPRester is None:  # pragma: no cover
                raise ImportError("pymatgen legacy MPRester is not available.")
            with MPRester(api_key) as m:
                structs = m.get_structures(chemsys_formula=formula)

        if len(structs) == 0:
            raise ValueError("Could not find composition in Materials Project Database, please supply a structure.")

        # Default to first found structure
        struct = structs[0]["structure"] if isinstance(structs[0], dict) else structs[0]  # type: ignore[index]

        if 0 not in (spec[1] for spec in sanit_species):  # If everything's charged
            if determine_oxi == "BV":
                bva = BVAnalyzer()
                struct = bva.get_oxi_state_decorated_structure(struct)

            elif determine_oxi == "comp_ICSD":
                comp = struct.composition
                oxi_transform = OxidationStateDecorationTransformation(comp.oxi_state_guesses(max_sites=-50)[0])  # type: ignore[union-attr]
                struct = oxi_transform.apply_transformation(struct)
                logger.info("Charge assigned based on ICSD statistics")

            elif determine_oxi == "both":
                try:
                    bva = BVAnalyzer()
                    struct = bva.get_oxi_state_decorated_structure(struct)
                    logger.info("Oxidation states assigned using bond valence")
                except (ValueError, RuntimeError):
                    comp = struct.composition
                    oxi_transform = OxidationStateDecorationTransformation(comp.oxi_state_guesses(max_sites=-50)[0])  # type: ignore[union-attr]
                    struct = oxi_transform.apply_transformation(struct)
                    logger.info("Oxidation states assigned based on ICSD statistics")
            else:
                raise ValueError(
                    f"Argument for 'determine_oxi', <{determine_oxi}> is not valid. Choose either 'BV','comp_ICSD' or 'both'."
                )
        lattice_mat = struct.lattice.matrix  # type: ignore[union-attr]

        lattice_param = 1.0  # Scaling factor; lattice_mat already contains actual vectors

        sites, _ = SmactStructure.__parse_py_sites(struct)

        return SmactStructure(
            _cast("list[tuple[str, int, int] | tuple[smact.Species, int]]", sanit_species),
            lattice_mat,
            sites,
            lattice_param,
            sanitise_species=False,
        )

    @staticmethod
    def from_file(fname: str) -> SmactStructure:
        """
        Create SmactStructure from a POSCAR file.

        Args:
        ----
            fname: The name of the POSCAR file.
                See :meth:`~.as_poscar` for format specification.

        Returns:
        -------
            :class:`~.SmactStructure`

        """
        with open(fname) as f:
            return SmactStructure.from_poscar(f.read())

    @staticmethod
    def from_poscar(poscar: str) -> SmactStructure:
        """
        Create SmactStructure from a POSCAR string.

        Args:
        ----
            poscar: A SMACT-formatted POSCAR string.
                See :meth:`~.as_poscar` for format specification.

        Returns:
        -------
            :class:`~.SmactStructure`

        """
        lines = poscar.split("\n")

        # Find stoichiometry
        total_specs = [int(x) for x in lines[6].split(" ")]
        hcf = reduce(gcd, total_specs)
        total_specs = [int(x / hcf) for x in total_specs]

        species = []
        for spec_str, stoic in zip(lines[0].split(" "), total_specs, strict=True):
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

        return SmactStructure(
            _cast("list[tuple[str, int, int] | tuple[smact.Species, int]]", species),
            lattice,
            sites,
            lattice_param,
        )

    def _format_style(
        self,
        template: str,
        delim: str | None = " ",
        include_ground: bool | None = False,
    ) -> str:
        """
        Format a given template string with the composition.

        Formats a python template string with species information,
        with each species separated by a given delimiter.

        Args:
        ----
            template: Template string to format, using python's
                curly brackets notation. Supported keywords are
                `ele` for the elemental symbol, `stoic` for the
                stoichiometry, `charge` for the absolute value
                of oxidation state and `sign` for the
                oxidation state's sign.
            delim: The delimiter between species' templates.
            include_ground: Whether to include the charge and sign
                of neutral species.

        Returns:
        -------
            String of templates formatted for each species, separated
                by `delim`.

        Examples:
        --------
            >>> s = SmactStructure.from_file("tests/files/CaTiO3.txt")
            >>> template = "{stoic}x{ele}{charge}{sign}"
            >>> print(s._format_style(template))
            1xCa2+ 3xO2- 1xTi4+

        """
        if include_ground:
            return delim.join(  # type: ignore[union-attr]
                template.format(
                    ele=specie[0],
                    stoic=specie[2],
                    charge=abs(specie[1]),
                    sign="+" if specie[1] >= 0 else "-",
                )
                for specie in self.species
            )

        return delim.join(  # type: ignore[union-attr]
            template.format(
                ele=specie[0],
                stoic=specie[2],
                charge=abs(specie[1]) if specie[1] != 0 else "",
                sign=get_sign(specie[1]),
            )
            for specie in self.species
        )

    @staticmethod
    def _get_ele_stoics(species: list[tuple[str, int, int]]) -> dict[str, int]:
        """
        Get the number of each element type in the compound, irrespective of oxidation state.

        Args:
        ----
            species: See :meth:`~.__init__`.

        Returns:
        -------
            eles: Dictionary of {element: stoichiometry}.

        Examples:
        --------
            >>> species = [("Fe", 2, 1), ("Fe", 3, 2), ("O", -2, 4)]
            >>> print(SmactStructure._get_ele_stoics(species))
            {'Fe': 3, 'O': 4}

        """
        eles = defaultdict(int)
        for specie in species:
            eles[specie[0]] += specie[2]

        return dict(eles)

    def has_species(self, species: tuple[str, int]) -> bool:
        """Determine whether a given species is in the structure."""
        return species in map(itemgetter(0, 1), self.species)

    def get_spec_strs(self) -> list[str]:
        """
        Get string representations of the constituent species.

        Returns:
        -------
            A list of strings, formatted as '{element}{charge}{sign}'.

        Examples:
        --------
            >>> s = SmactStructure.from_file("tests/files/CaTiO3.txt")
            >>> s.get_spec_strs()
            ['Ca2+', 'O2-', 'Ti4+']

        """
        return self._format_style("{ele}{charge}{sign}").split(" ")

    def composition(self) -> str:
        """
        Generate a key that describes the composition.

        Key format is '{element}_{stoichiometry}_{charge}{sign}'
        with no delimiter, *sans brackets*. Species are ordered as stored within
        the structure, see :class:`~.SmactStructure`.

        Returns:
        -------
            Key describing constituent species.

        Examples:
        --------
            >>> s = SmactStructure.from_file("tests/files/CaTiO3.txt")
            >>> print(s.composition())
            Ca_1_2+O_3_2-Ti_1_4+

        """
        comp_style = "{ele}_{stoic}_{charge}{sign}"
        return self._format_style(comp_style, delim="", include_ground=True)

    def as_poscar(self) -> str:
        """
        Represent the structure as a POSCAR file compatible with VASP5.

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
        -------
            str: POSCAR-style representation of the structure.

        """
        poscar = " ".join(self.get_spec_strs()) + "\n"

        poscar += f"{self.lattice_param}\n"

        poscar += "\n".join(" ".join(map(str, vec)) for vec in self.lattice_mat.tolist()) + "\n"

        spec_count = {spec: len(coords) for spec, coords in self.sites.items()}

        poscar += self._format_style("{ele}") + "\n"

        poscar += " ".join(str(spec_count[spec]) for spec in self.get_spec_strs()) + "\n"

        poscar += "Cartesian\n"
        for spec, coords in self.sites.items():
            for coord in coords:
                poscar += " ".join(map(str, coord))
                poscar += f" {spec}\n"

        return poscar

    def as_py_struct(self) -> pmg_Structure:
        """
        Represent the structure as a pymatgen Structure object.

        Returns:
        -------
            pmg_Structure: pymatgen Structure object.

        """
        return pmg_Structure.from_str(self.as_poscar(), fmt="poscar")

    def reduced_formula(self) -> str:
        """
        Generate a reduced formula for the structure.

        Returns:
        -------
            str: Reduced formula of the structure.

        Examples:
        --------
            >>> s = SmactStructure.from_file("tests/files/CaTiO3.txt")
            >>> print(s.reduced_formula())
            CaTiO3

        """
        return self.as_py_struct().composition.reduced_formula
