#!/usr/bin/env python 
# Using the ase spacegroup module this can build the structure, from 
# the composition, as defined in the smact_lattice module. 
#TODO:
# Estimate the cell parameters based on radii from tables.
# Add further types, Spinel, Flourite, Delafossite ....

from ase.lattice.spacegroup import crystal


def cubic_perovskite(species):
	 system = crystal((species), 
	 [(0,0,0), (0.5, 0.5, 0.5), (0.5, 0.5, 0)],
         spacegroup=221, cellpar=[6,6,6,90,90,90])
	
         return system
