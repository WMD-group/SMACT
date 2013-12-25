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
import ase
import smact_builder as builder
from pyspglib import spglib
from ase.lattice.spacegroup import Spacegroup
import numpy as np

#------------------------------------------------------------------------------------------
def get_sg(lattice):
    # Get the space-group of the system
    spacegroup = spglib.get_spacegroup(lattice, symprec=1e-5)
    space_split=spacegroup.split()
    spg_num = space_split[1].replace('(','').replace(')','')
    sg = Spacegroup(int(spg_num))
    return sg
#------------------------------------------------------------------------------------------
def are_eq(A,B,tolerance=1e-4):
    """Check two arrays for tolearnce [1,2,3]==[1,2,3]; but [1,3,2]!=[1,2,3]"""
    are_eq = True
    if len(A) != len(B):
	are_eq = False
    else:
        i = 0
        while i < len(A):
	    if abs(A[i] - B[i]) > tolerance:
	        are_eq = False
            i = i + 1
    return are_eq
#------------------------------------------------------------------------------------------
def get_inequivalent_sites(sub_lattice, lattice):
    """Given a sub lattice, returns symmetry unique sites for substitutions"""
    sg = get_sg(lattice)
    inequivalent_sites = []
    i = 0
    for site in sub_lattice:
        new_site = True
# Check against the existing members of the list of inequivalent sites
        for inequiv_site in inequivalent_sites:
	    if are_eq(site, inequiv_site) == True:
	        new_site = False
# Check against symmetry related memebers of the list of inequivalent sites
            equiv_inequiv_sites,junk = sg.equivalent_sites(inequiv_site)    	
	    for equiv_inequiv_site in equiv_inequiv_sites:
	        if are_eq(site, equiv_inequiv_site) == True:
	   	    new_site = False
        if new_site == True:
	    inequivalent_sites.append(site)
        i = i + 1
    return inequivalent_sites
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
