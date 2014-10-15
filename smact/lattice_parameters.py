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

def wurtzite(shannon_radius):
    '''The lattice parameters of the wurtzite structure
    Args:
	shannon_radius : a list containing the radii of the a,b ions
    Returns:
	a,b,c : real number values of the lattice constants
	alpha,beta,gamma : real number values of the lattice angles
    '''
    shannon_radius.sort(reverse=True) # Geometry assumes atom A is larger
    # "Ideal" wurtzite structure
    # c/a = 1.633, u = 0.375
    #
    alpha = 90
    beta = 90
    gamma = 120
    # 
    # Scenario A: A atoms are touching
    #   i.e. height is that of two tetrahegons with side length a
    #    = 2 * sqrt(2/3) * a
    if shannon_radius[0] > 0.817*(shannon_radius[0]+shannon_radius[1]):
        a = 2 * shannon_radius[0]
        b = a
        c = 2 * np.sqrt(2./3.) * a
    else:
        # Scenario B: regular wurtzite, similar sizes
        a = 2*0.817*(shannon_radius[0]+shannon_radius[1])  # 0.817 is sin(109.6/2)
        b = a
        c = (shannon_radius[0]+shannon_radius[1])*(2+2*0.335)  # 0.335 is sin(109.6-90)
#    inner_space = a * (6**0.5) - (4*shannon_radius[0])
    return a,b,c,alpha,beta,gamma
