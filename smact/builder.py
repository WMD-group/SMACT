"""A collection of functions for building certain lattice types."""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

from ase.spacegroup import crystal

if TYPE_CHECKING:
    from ase import Atoms

from smact.lattice import Lattice, Site


def _tile_oxidation_states(
    base_oxidation_states: list[list[int]],
    n_sites: int,
) -> list[list[int]]:
    """Tile the base unit-cell oxidation states to match the number of sites in an expanded cell."""
    n_base = len(base_oxidation_states)
    if n_sites == n_base:
        return base_oxidation_states
    if n_sites % n_base != 0:
        msg = (
            f"Number of sites ({n_sites}) is not a multiple of the base "
            f"oxidation states pattern ({n_base})"
        )
        raise ValueError(msg)
    return list(itertools.islice(itertools.cycle(base_oxidation_states), n_sites))


def cubic_perovskite(
    species: list[str],
    cell_par: list[float] | None = None,
    repetitions: list[int] | None = None,
) -> tuple[Lattice, Atoms]:
    """
    Build a perovskite cell using the crystal function in ASE.

    Args:
    ----
        species (str): Element symbols
        cell_par (list): Six floats/ints specifying 3 unit cell lengths and 3 unit cell angles.
        repetitions (list): Three floats specifying the expansion of the cell in x,y,z directions.

    Returns:
    -------
        SMACT Lattice object of the unit cell,
        ASE crystal system of the unit cell.

    """
    if repetitions is None:
        repetitions = [1, 1, 1]
    if cell_par is None:
        cell_par = [6, 6, 6, 90, 90, 90]
    system = crystal(
        (species),
        basis=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0)],
        spacegroup=221,
        size=repetitions,
        cellpar=cell_par,
    )

    base_oxidation_states = [[2]] + [[4]] + [[-2]] * 3
    oxidation_states = _tile_oxidation_states(base_oxidation_states, len(system))
    sites_list = [Site(pos, ox) for pos, ox in zip(system.get_scaled_positions(), oxidation_states, strict=True)]

    return Lattice(sites_list), system


def wurtzite(
    species: list[str],
    cell_par: list[float] | None = None,
    repetitions: list[int] | None = None,
) -> tuple[Lattice, Atoms]:
    """
    Build a wurzite cell using the crystal function in ASE.

    Args:
    ----
        species (str): Element symbols
        cell_par (list): Six floats/ints specifying 3 unit cell lengths and 3 unit cell angles.
        repetitions (list): Three floats specifying the expansion of the cell in x,y,z directions.

    Returns:
    -------
        SMACT Lattice object of the unit cell,
        ASE crystal system of the unit cell.

    """
    if repetitions is None:
        repetitions = [1, 1, 1]
    if cell_par is None:
        cell_par = [3, 3, 6, 90, 90, 120]
    system = crystal(
        (species),
        basis=[(2.0 / 3.0, 1.0 / 3.0, 0), (2.0 / 3.0, 1.0 / 3.0, 5.0 / 8.0)],
        spacegroup=186,
        size=repetitions,
        cellpar=cell_par,
    )

    base_oxidation_states = [[2]] * 2 + [[-2]] * 2
    oxidation_states = _tile_oxidation_states(base_oxidation_states, len(system))
    sites_list = [Site(pos, ox) for pos, ox in zip(system.get_scaled_positions(), oxidation_states, strict=True)]
    return Lattice(sites_list), system
