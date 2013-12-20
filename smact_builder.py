#!/usr/bin/env python 
# Using the ase spacegroup module this can build the structure, from 
# the composition, as defined in the smact_lattice module. 
#TODO:
# Estimate the cell parameters based on radii from tables.
# Add further types, Spinel, Flourite, Delafossite ....

################################################################################
# Copyright Keith Butler (2013)                                                #
#                                                                              #
#  This file is part of SMACT: smact_builder.py is free software: you can      #
#  redistribute it and/or modify it under the terms of the GNU General Public  #
#  License as published by the Free Software Foundation, either version 3 of   #
#  the License, or (at your option) any later version.                         #
#  This program is distributed in the hope that it will be useful, but WITHOUT #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   #
#  more details.                                                               #
#  You should have received a copy of the GNU General Public License along with#
#  this program.  If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

from ase.lattice.spacegroup import crystal


def cubic_perovskite(species):
	 system = crystal((species), 
	 [(0,0,0), (0.5, 0.5, 0.5), (0.5, 0.5, 0)],
         spacegroup=221, cellpar=[6,6,6,90,90,90])
	
         return system
