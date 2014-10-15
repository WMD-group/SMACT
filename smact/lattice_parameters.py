#!/usr/bin/env python 
# This script can calculate roughly the lattice parameters of a lattice type, based 
# on the radii of the species on each site.
################################################################################
# Copyright Tim Gauntlett   (2014)                                                         #
#                                                                              #
# This file is part of SMACT: parameters.py is free software: you can          #
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

def perovskite_parameters(shannon_radius): #Cubic Pervoskite
    '''The lattice parameters of the cubic perovskite structure
    Args:
	shannon_radius : a list containing the radii of the a,b,c ions
    Returns:
	a,b,c : real number values of the lattice constants
	alpha,beta,gamma : real number values of the lattice angles
	'''
    limiting_factors=[2*sum(shannon_radius[1:])]
    a = max(limiting_factors)
    b = a
    c = a
    space = a * np.sqrt(3) - 2 * shannon_radius[1]
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,space,alpha,beta,gamma


