#!/usr/bin/env python

################################################################################
# Copyright Daniel Davies (2013)                                               #
#                                                                              #
#  This file is part of SMACT: ValenceMax.py is free software: you can         #
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

from numpy import product,sqrt
from smact_data import get_eig,get_covalent
from scipy.constants import hbar, m_e
from scipy.constants import physical_constants
J2eV = physical_constants['joule-electron volt relationship'][0]

"""Get element names"""
An = raw_input("Enter Anion Symbol: ")
Cat = raw_input("Enter Cation Symbol: ")

"""calculations"""
V2 = (2.16*(hbar)**2)/(
                       m_e* (
                             get_covalent(An)*1E-10 +
                             get_covalent(Cat)*1E-10
                             )**2
                       )

V3 =  (get_eig(Cat) - get_eig(An))/2

Ev = ((get_eig(Cat) + get_eig(An))/2) - (sqrt((V2**2)+(V3**2)))


print An, "Eigenvalue:", get_eig(An)
print Cat, "Eigenvalue:", get_eig(Cat)
print " "
print "V2:", V2*J2eV, "eV"
print "V3:", V3, "eV"
print "Valence Band Minimum = ", Ev, "eV"
