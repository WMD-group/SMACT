#!/usr/bin/env python
# Using the ase spacegroup module this can build the structure, from
# the composition, as defined in the smact_lattice module.
# TODO:
# Estimate the cell parameters based on radii from tables.
# Add further types, Spinnel, Flourite, Delafossite ....

# Implement Structure class, c.f. dev_docs.

from ase.spacegroup import crystal

from smact.lattice import Lattice, Site


def cubic_perovskite(
    species, cell_par=[6, 6, 6, 90, 90, 90], repetitions=[1, 1, 1]
):
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
    system = crystal(
        (species),
        basis=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0)],
        spacegroup=221,
        size=repetitions,
        cellpar=cell_par,
    )

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
    system = crystal(
        (species),
        basis=[(2.0 / 3.0, 1.0 / 3.0, 0), (2.0 / 3.0, 1.0 / 3.0, 5.0 / 8.0)],
        spacegroup=186,
        size=repetitions,
        cellpar=[3, 3, 6, 90, 90, 120],
    )

    sites_list = []
    oxidation_states = [[1], [2], [3], [4]] + [[-1], [-2], [-3], [-4]]

    for site in zip(system.get_scaled_positions(), oxidation_states):
        sites_list.append(Site(site[0], site[1]))
    return Lattice(sites_list, oxidation_states), system
