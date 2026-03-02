"""
smact.distorter: Module for generating symmetry-unique substitutions on a given sub-lattice.

As input it takes the ASE crystal object (as built by smact.builder)
and the sub-lattice on which substitutions are to be made.
"""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING, cast

import smact

if TYPE_CHECKING:
    from ase import Atoms
    from ase.atom import Atom

try:
    from pyspglib import spglib  # type: ignore[import-not-found]
except ImportError:
    try:
        import spglib
    except ImportError as e:  # pragma: no cover
        msg = "Could not load spglib. Install it with: pip install spglib"
        raise ImportError(msg) from e

from ase.spacegroup import Spacegroup


def get_sg(lattice: Atoms) -> Spacegroup:
    """
    Get the space-group of the system.

    Args:
    ----
        lattice: the ASE crystal class
    Returns:
        sg (int): integer number of the spacegroup

    """
    spacegroup: str | None = spglib.get_spacegroup(lattice, symprec=1e-5)  # type: ignore[arg-type]  # spglib expects Cell, ASE Atoms is compatible at runtime
    if spacegroup is None:
        msg = "spglib could not determine the spacegroup of the lattice"
        raise ValueError(msg)
    space_split = spacegroup.split()
    spg_num = space_split[1].replace("(", "").replace(")", "")
    return Spacegroup(int(spg_num))


def get_inequivalent_sites(
    sub_lattice: list[list[float]],
    lattice: Atoms,
) -> list[list[float]]:
    """
    Given a sub lattice, returns symmetry unique sites for substitutions.

    Args:
    ----
        sub_lattice (list of lists): array containing Cartesian coordinates
            of the sub-lattice of interest

        lattice (ASE crystal): the total lattice

    Returns:
    -------
        List of sites

    """
    sg = get_sg(lattice)
    inequivalent_sites = []
    for site in sub_lattice:
        new_site = True
        # Check against the existing members of the list of inequivalent sites
        if len(inequivalent_sites) > 0:
            for inequiv_site in inequivalent_sites:
                if smact.are_eq(site, inequiv_site):
                    new_site = False
                # Check against symmetry related members of the list of inequivalent sites
                equiv_inequiv_sites, _ = sg.equivalent_sites(inequiv_site)
                for equiv_inequiv_site in equiv_inequiv_sites:
                    if smact.are_eq(site, equiv_inequiv_site):
                        new_site = False

        if new_site:
            inequivalent_sites.append(site)

    return inequivalent_sites


def make_substitution(lattice: Atoms, site: list[float], new_species: str) -> Atoms:
    """
    Change atomic species on lattice site to new_species.

    Args:
    ----
        lattice (ASE crystal): Input lattice
        site (list): Cartesian coordinates of the substitution site
        new_species (str): New species

    Returns:
    -------
        lattice

    """
    # deepcopy is necessary, otherwise changes applied to the clone also apply to the parent object.
    new_lattice = copy.deepcopy(lattice)
    lattice_sites = new_lattice.get_scaled_positions()
    for i, lattice_site in enumerate(lattice_sites):
        if smact.are_eq(lattice_site, site):
            atom = cast("Atom", new_lattice[i])
            atom.symbol = new_species
    return new_lattice


def build_sub_lattice(lattice: Atoms, symbol: str) -> list[list[float]]:
    """
    Generate a sub-lattice of the lattice based on equivalent atomic species.

    Args:
    ----
        lattice (ASE crystal class): Input lattice
        symbol (string): Symbol of species identifying sub-lattice

    Returns:
    -------
        list of lists:
            sub_lattice: Cartesian coordinates of the sub-lattice of symbol

    """
    sub_lattice = []
    atomic_labels = lattice.get_chemical_symbols()
    positions = lattice.get_scaled_positions()
    for atom, pos in zip(atomic_labels, positions, strict=True):
        if atom == symbol:
            sub_lattice.append(pos)
    return sub_lattice
