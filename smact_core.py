#!/usr/bin/env python

################################################################################
#  Copyright Adam J. Jackson, Daniel Davies (2013)                             #
#                                                                              #
#  This file is part of SMACT: smact_core.py is free software: you can         #
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

"""
Core classes and functions for SMACT
"""

class element(object):
    """Class providing standard chemical data for elements."""
    def __init__(self, symbol):
                """
                Collection of standard elemental properties for given element.
                Data is drawn from "data/element.txt", part of the Open Babel package.
                element = element('Symbol')

                Methods
                -------
                element.symbol: Elemental symbol used to retrieve data

                element.name: Full name of element

                element.covalent_radius: Covalent radius in AA (1.6 if unknown)
                
                element.vdw_radius: van der Waals radius in AA (2.0 if unknown)

                element.pauling_eneg: Pauling electronegativity (0.0 if unknown)

                element.ionpot: Ionisation potential in eV (0.0 if unknown)

                element.e_affinity: Eletron affinity in eV (0.0 if unknown)                
                """
                with open('data/element.txt','r') as f:
                        data = f.readlines()
                elementdata = 0
                for line in data:
                        if not line.startswith('#'):
                                l = line.split()
                                if (l[1] == symbol):
                                        elementdata = l
                                        break
                if not elementdata:
                        raise NameError('Element {0} not found'.format(symbol))
                else:
                        self.symbol=          str(symbol)
                        self.name=            str(elementdata[14])
                        self.covalent_radius= float(elementdata[3])
                        self.vdw_radius=      float(elementdata[5])
                        self.pauling_eneg=    float(elementdata[8])
                        self.ionpot=          float(elementdata[9])
                        self.e_affinity =     float(elementdata[10])
#                       self.eigenval = 

        
