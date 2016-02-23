#!/usr/bin/env python

################################################################################
# Copyright Daniel Davies (2013)                                               #
#                                                                              #
# This file is part of SMACT: Band_gap_simple.py is free software: you can     #
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
import smact
from numpy import sqrt


def band_gap_simple(verbose=False, anion=None, cation=None,
                    distance=None, elements_dict=None):

    """
    Estimates the band gap from elemental data.

    The band gap is estimated using the principles outlined in
    Harrison's 1980 work "Electronic Structure and the Properties of
    Solids: The Physics of the Chemical Bond".

    Args:
        Anion (str): Element symbol of the dominant anion in the system.

        Cation (str): Element symbol of the the dominant cation in the system.
        Distance (float or str): Nuclear separation between anion and cation
                i.e. sum of ionic radii.
        verbose: An optional True/False flag. If True, additional information is
                printed to the standard output. [Defult: False]
        elements_dict (dict): Element symbol keys with smact.Element
                objects values. This may be provided to prevent
                excessive reading of data files.

    Returns (float): Band gap in eV.

    """

    # Set constants
    hbarsq_over_m = 7.62

    # Get anion and cation
    if anion:
        An = anion
    if cation:
        Cat = cation
    if not anion:
        An = raw_input("Enter Anion Symbol:")
    if not cation:
        Cat = raw_input("Enter Cation Symbol:")
    if distance:
        d = float(distance)
    if not distance:
        d = float(raw_input("Enter internuclear separation (Angstroms): "))

    # Get elemental data:
    if not elements_dict:
        elements_dict = smact.element_dictionary((An, Cat))
    An, Cat = elements_dict[An], elements_dict[Cat]

    # Calculate values of equation components
    V1_Cat = (Cat.eig - Cat.eig_s)/4
    V1_An = (An.eig - An.eig_s)/4
    V1_bar = (V1_An + V1_Cat)/2
    V2 = 2.16 * hbarsq_over_m / (d**2)
    V3 = (Cat.eig - An.eig)/2
    alpha_m = (1.11*V1_bar)/sqrt(V2**2 + V3**2)

    # Calculate Band gap [(3-43) Harrison 1980 ]
    Band_gap = (3.60/3)*(sqrt(V2**2 + V3**2))*(1-alpha_m)
    if verbose:
        print "V1_bar = ", V1_bar
        print "V2 = ", V2
        print "alpha_m = ", alpha_m
        print "V3 = ", V3

    return Band_gap

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description ="""Compound band gap estimates from elemental data.""")
    parser.add_argument("-a", "--anion", type=str,
                        help ="""Element symbol for anion.""")
    parser.add_argument("-c", "--cation", type=str,
                        help ="""Element symbol for cation.""")
    parser.add_argument("-d", "--distance", type=float,
                        help ="""Internuclear separation.""")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help ="""More Verbose output.""")
    args = parser.parse_args()

    print band_gap_simple(verbose=args.verbose,anion=args.anion,
                          cation=args.cation,distance=args.distance)
