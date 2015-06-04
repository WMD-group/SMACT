#!/usr/bin/env python

################################################################################
# Copyright Daniel Davies, Adam J. Jackson (2013)                              #
#                                                                              #
# This file is part of SMACT: compound_electroneg.py is free software: you can #
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

import sys
from numpy import product
from smact.data import get_mulliken

def compound_electroneg(verbose=False,elements=None,stoichs=None):

    """Estimate Mulliken electronegativity of compound from elemental data.
		Uses get_mulliken function which uses elemental ionisation potentials
		and electron affinities. 

    Geometric mean is used (n-th root of product of components), e.g.:

    X_Cu2S = (X_Cu * X_Cu * C_S)^(1/3)

    Args:
        elements: A list of elements given as standard elemental symbols.
		Optional: if not used, interactive input of space separated
		elemental symbols will be offered.
        stoichs: A list of stoichiometries, given as integers or floats.
		Optional: if not used, interactive input of space separated
		integers  will be offered.
        verbose: An optional True/False flag. If True, additional information is
               printed to the standard output. [Default: False]

    Returns:
        Electronegativity: Estimated electronegativity as a float.
                         Electronegativity is a dimensionless property.

    Raises:
        (There are no special error messages for this function.)
    
    """
    if elements:
	elementlist = elements
    if stoichs:
	stoichslist = stoichs
    
    """Get elements and stoichiometries if not provided as argument"""
    if not elements:
        elementlist = list(raw_input("Enter elements (space separated): ").split(" "))	    
    if not stoichs:
        stoichslist = list(raw_input("Enter stoichiometries (space separated): ").split(" "))

  
    """Convert stoichslist from string to float"""
    stoichslist=map(float, stoichslist)

    """Get mulliken values for each element"""
    for i in range(0,len(elementlist)):
        elementlist[i]=get_mulliken(elementlist[i])

    """Print optional list of element electronegativities.
    This may be a useful sanity check in case of a suspicious result.
    """
    if verbose:
        print "Electronegativities of elements=", elementlist

    """Raise each electronegativity to its appropriate power
    to account for stoichiometry.
    """
    for i in range(0,len(elementlist)):
        elementlist[i]=[elementlist[i]**stoichslist[i]]
    """Calculate geometric mean (n-th root of product)"""
    prod = product(elementlist)
    compelectroneg = (prod)**(1.0/(sum(stoichslist)))

    """Print optional formatted output."""
    if verbose:
        print "Geometric mean = Compound 'electronegativity'=", compelectroneg
      
    return compelectroneg

"""Wrapper for command line usage: argparse passes command line arguments
to python and proves help function.
"""

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Compound electronegativity from geometric mean of
                       elemental Mulliken electronegativities.""")
    parser.add_argument("-e", "--elements", type=str,
                        help="""Space-separated string of elements
                                (e.g. \"Cu Zn Sn S\")"""
                        )
    parser.add_argument("-s", "--stoichiometry", type=str,
                        help="""Space-separated string of stoichiometry
                                (e.g. \"2 1 1 4\")""" % ()
                        )
    parser.add_argument("-v", "--verbose", help="More verbose output [default]",
                        action="store_true")
    parser.add_argument("-q", "--quiet", help="Quiet output",
                       action="store_true")
    args=parser.parse_args()
    if args.quiet:
        verbose_flag=False
    else:
        verbose_flag=True

    print compound_electroneg(verbose=verbose_flag,elements=args.elements,
                              stoichs=args.stoichiometry)

