#!/usr/bin/env python
# Using the ase spacegroup module this can build the structure, from
# the composition, as defined in the smact_lattice module.
#TO DO:
# Calculating the lattice parameters from the elemnts involved
################################################################################
# Copyright Tim Gauntlett   (2014)  				                           #
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

from __future__ import print_function
from __future__ import division
from builtins import input
from past.utils import old_div
import smact
import numpy as np
import sys
import csv

def main():
    crystal_elements = input('Which elements? Separate chemical symbols with a space. ')
    crystal_elements = crystal_elements.split()

    # Processing Lattices
    #------------------------------------------------------------------------------------------

    ## A-types
    print (" \n=========================================")
    print ("    AAA         ####### ##   ## #####  #####")
    print ("   AA A           ##     ## ## ##  ## ##")
    print ("  AAAAAA  ####   ##      ##   #####  #####")
    print (" AA   AA        ##      ##   ##     ##")
    print ("AA     AA      ##      ##   ##     #####")
    print ("=========================================")
    print (" \n ### A type lattices ###\n ")

    for lattice in a_lattices:
        print((' \n===============\n ') + lattice.__name__.upper() + (' Lattice \n===============\n  '))
        for elements in crystal_elements:
            x = smact.Element(elements)
            a,b,c,alpha,beta,gamma = lattice(x.covalent_radius)
            print(elements,"%.2f"%a,"%.2f"%b,"%.2f"%c,"%.1f"%alpha,"%.1f"%beta,"%.1f"%gamma)
    print(' \n Key :: Element a b c alpha beta gamma \n///////////////////////////////////////////////////////////////////////////////////////\n ')

    oxidation =[]
    coordination =[]

    i=0
    for elements in crystal_elements:
        poss_ox = []
        poss_cood = []
        with open('data/shannon_radii.csv','rU') as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == elements:
                    poss_ox.append(row[1])
                    poss_cood.append(row[2])
        poss_ox_clean = [] ### Duplicates removed
        for ox in poss_ox:
            if ox not in poss_ox_clean:
                poss_ox_clean.append(ox)
        print(elements)
        print(poss_ox_clean)
        if len(poss_ox_clean) == 1:
            oxidation.append(int(poss_ox_clean[0]))
        else:
            oxidation.append(int(input('Oxidation state of ' + elements + ' ')))
        poss_co_clean = [] ### Duplicates removed
        with open('data/shannon_radii.csv','rU') as f:
            reader = csv.reader(f)
            poss_cood = []
            for row in reader:
                if row[0] == elements and int(row[1]) == int(oxidation[-1]):
                    poss_cood.append(row[2])
        for co in poss_cood:
            if co not in poss_co_clean:
                poss_co_clean.append(co)
        print(elements)
        print(poss_co_clean)
        if len(poss_co_clean) == 1:
            coordination.append(poss_co_clean[0])
        else:
            coordination.append(input('Coordination state of ' + elements + ' '))

        print(elements,oxidation[i],coordination[i])
        i+=1

    ## B-types
    #-------------------------------------------------------------------------------------#
    print (" \n=========================================")
    print ("    BBBB     ####### ##   ## #####  #####")
    print ("   BB BB       ##     ## ## ##  ## ##")
    print ("  BBBB   ###  ##      ##   #####  #####")
    print (" BB BB       ##      ##   ##     ##")
    print ("BBBB        ##      ##   ##     #####")
    print ("=========================================")

    print (" \n ### B type lattices ###\n ")

    for lattice in b_lattices:
        print((' \n===============\n ') + lattice.__name__.upper() + (' Lattice \n===============\n  '))
        shannon_radius = []
        for i, elements in enumerate(crystal_elements):
            print(elements,oxidation[i],coordination[i])
            x = smact.Species(elements,oxidation[i],coordination[i])
            shannon_radius.append(float(x.shannon_radius))
        a,b,c,alpha,beta,gamma = lattice(shannon_radius)
        if lattice.__name__ == 'wurtzite':
            inner_space = a * (6**0.5) - (4*shannon_radius[0])
            print ("With a gap in the middle of " + "%.2f"%inner_space+" (diameter)")
        print(crystal_elements[0:2],"%.2f"%a,"%.2f"%b,"%.2f"%c,"%.1f"%alpha,"%.1f"%beta,"%.1f"%gamma)

    # E Types
    #-------------------------------------------------------------------------------------#
    print (" \n=========================================")
    print ("    EEEEE   ####### ##   ## #####  #####")
    print ("   EE         ##     ## ## ##  ## ##")
    print ("  EEEE  ###  ##      ##   #####  #####")
    print (" EE         ##      ##   ##     ##")
    print ("EEEEE      ##      ##   ##     #####")
    print ("=========================================")

    print (" \n ### E type lattices ###\n ")

    for lattice in e_lattices:
        print((' \n===============\n ') + lattice.__name__.upper() + (' Lattice \n===============\n  '))
        shannon_radius = []
        for i, elements in enumerate(crystal_elements):
            print(elements,oxidation[i],coordination[i])
            x = smact.Species(elements,oxidation[i],coordination[i])
            shannon_radius.append(float(x.shannon_radius))
        a,b,c,space,alpha,beta,gamma = lattice(shannon_radius)
        print(crystal_elements[0:],"%.2f"%a,"%.2f"%b,"%.2f"%c,"%.1f"%alpha,"%.1f"%beta,"%.1f"%gamma)
        print ("With a gap in the middle of " + "%.2f"%space+" (diameter)")

if __name__ == '__main__':
    main()
