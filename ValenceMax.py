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
"get_eig function not yet working"""
from smact_data import get_eig,get_covalent
from scipy.constants import hbar m_e

"""Get element names"""
An = 
Cat =

"""calculations"""
V2 = (2.16*(hbar)**2)/m_e*((get_covalent(An))+(get_covalent(Cat))**2)

V3 =  (get_eig(Cat) + get_eig(An))/2 

Ev = V3 + sqrt((V2**2)+(V3**2))

print "Valence Band Minimum = ", Ev
