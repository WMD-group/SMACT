################################################################################
# Copyright Keith T Butler    (2013)                                           #
#                                                                              #
# This file is part of SMACT: smact_distorter.py is free software: you can     #
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

"""
smact_distorter: Module for generating symmetry-unique substitutions on a given sub-lattice.

As input it takes the ASE crystal object (as built by smact_builder)
and the sub-lattice on which substitutions are to be made. 
There is an example of how to use the code in Example_distort.py

---------------------------------------------------------------------
TODO:
Add a functionality to check two Atoms objects against one another 
for equivalence.
---------------------------------------------------------------------

"""

import ase
import copy
import smact_builder as builder
import smact_core as core
from pyspglib import spglib
from ase.lattice.spacegroup import Spacegroup
import numpy as np

#------------------------------------------------------------------------------------------
def get_sg(lattice):

    """
    Get the space-group of the system

    Args:
	lattice: the ASE crystal class
	Returns:
	sg: integer number of the spacegroup

    """
    spacegroup = spglib.get_spacegroup(lattice, symprec=1e-5)
    space_split=spacegroup.split()
    spg_num = space_split[1].replace('(','').replace(')','')
    sg = Spacegroup(int(spg_num))
    return sg
#------------------------------------------------------------------------------------------
def get_inequivalent_sites(sub_lattice, lattice):
    """Given a sub lattice, returns symmetry unique sites for substitutions

	Args:
	sub_lattice: array containing Cartesian coordinates of the sub-lattice of interest
	lattice: ASE crystal class, the total lattice
     """
    sg = get_sg(lattice)
    inequivalent_sites = []
    i = 0
    for site in sub_lattice:
        new_site = True
# Check against the existing members of the list of inequivalent sites
        for inequiv_site in inequivalent_sites:
	    if core.are_eq(site, inequiv_site) == True:
	        new_site = False
# Check against symmetry related members of the list of inequivalent sites
            equiv_inequiv_sites,junk = sg.equivalent_sites(inequiv_site)    	
	    for equiv_inequiv_site in equiv_inequiv_sites:
	        if core.are_eq(site, equiv_inequiv_site) == True:
	   	    new_site = False
        if new_site == True:
	    inequivalent_sites.append(site)
        i = i + 1
    return inequivalent_sites
#------------------------------------------------------------------------------------------
def make_substitution(lattice,site,new_species):
    """Change the atomic species @ site in lattice to new_species [atomic number]

    Args:
	lattice: ASE crystal class
	site: array, containing the Cartesian coordinates of the substitution site
	new_species: string, the new species

     """
    i = 0
# NBNBNBNB  It is necessary to use deepcopy for objects, otherwise changes applied to a clone
# will also apply to the parent object.
    new_lattice = copy.deepcopy(lattice)
    lattice_sites = new_lattice.get_scaled_positions()
    for lattice_site in lattice_sites:
	if core.are_eq(lattice_site, site):
	    new_lattice[i].symbol = new_species
	i = i + 1
    return new_lattice
#------------------------------------------------------------------------------------------
def build_sub_lattice(lattice,symbol):
    """Generate a sub-lattice of the lattice based on equivalent atomic species

    Args:
        lattice: ASE crystal class
	symbol: The name of the species whose sub-lattice you wish to build

    Returns:
	sub_lattice: Cartesian coordinates of the sub-lattice of symbol

    """

    sub_lattice = []
    i = 0
    atomic_labels = lattice.get_chemical_symbols()    
    positions = lattice.get_scaled_positions()
    for atom in atomic_labels:
        if atom == symbol:
            sub_lattice.append(positions[i])
        i = i + 1
    return sub_lattice
#------------------------------------------------------------------------------------------


