#!/usr/bin/env python

################################################################################
# Copyright Daniel Davies, Adam J. Jackson (2013)                              #
#                                                                              #
#  This file is part of SMACT: data.py is free software: you can               #
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

import smact

def get_mulliken(element):
    """Get Mulliken electronegativity from the IE and EA

    Arguments:
        symbol (smact.Element or str): Element object or symbol

    Returns:
        mulliken (float): Mulliken electronegativity

    """
    if type(element) == str:       
        element = smact.Element(element)
    elif type(element) != smact.Element:
        raise Exception("Unexpected type: {0}".format(type(element)))        

    mulliken = (element.ionpot+element.e_affinity)/2.0

    return mulliken


################################################################################
#### The following functions are deprecated. New code should access these   ####
#### properties directly by forming an "element" object and calling the     ####
#### desired property, as in the wrappers below.                            ####
################################################################################

def get_pauling(element):
    """Pauling electronegativity of specified element.

    Arguments:
        symbol (smact.Element or str): Element object or symbol

    Returns:
        pauling_eneg (float): Pauling electronegativity

    """
    if type(element) == str:       
        element = smact.Element(element)
    elif type(element) != smact.Element:
        raise Exception("Unexpected type: {0}".format(type(element)))        

    return element.pauling_eneg

def get_covalent(symbol):
    """Covalent radius of specified element.
    Drawn from Open Babel data table.

    Arguments:
        symbol (string): Element label

    Returns:
        covalent_radius (float): Covalent radius
    """
    from smact import Element
    A = Element(symbol)
    return A.covalent_radius

def get_eig(symbol):
    """Eigenvalue of specified element.

    Arguments:
        symbol (string): Element label

    Returns:
        eig (float): Eigenvalue
    """    
    from smact import Element
    A=Element(symbol)
    return A.eig

def get_eig_s(symbol):
    """Eigenvalue of s-orbital in specified element.

    Arguments:
        symbol (string): Element label

    Returns:
        eig_s (float): Eigenvalue
    """        
    from smact import Element
    A=Element(symbol)
    return A.eig_s
