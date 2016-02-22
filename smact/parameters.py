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

import smact
import numpy as np

import sys
def query_yes_no(question, default="yes"):
        valid = {"yes": True, "y": True, "yer":True, "sure":True, "i am":True, "ye": True, "nope":False, "no": False, "n": False}
        if default is None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)

        while True:
            sys.stdout.write(question + prompt)
            choice = raw_input().lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                sys.stdout.write("Please respond with 'yes' or 'no' "
                                 "(or 'y' or 'n').\n")


# A-type lattices
#------------------------------------------------------------------------------------------

#A1#
def fcc(covalent_radius):
    a = 2 * 2**0.5 * covalent_radius
    b = 2 * 2**0.5 * covalent_radius
    c = 2 * 2**0.5 * covalent_radius
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma

#A2#
def bcc(covalent_radius):
    a = 4 * covalent_radius / np.sqrt(3)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma

#A3#
def hcp(covalent_radius):
    a = 2 * covalent_radius
    b = a
    c = (4./3.) * 6**0.5 * covalent_radius
    alpha = 90
    beta = 90
    gamma = 120
    return a,b,c,alpha,beta,gamma

#A4#
def diamond(covalent_radius):
    a = 8 * covalent_radius / np.sqrt(3)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma

#A5#
def bct(covalent_radius):
    a = 3.86 * covalent_radius
    b = a
    c = 2 * covalent_radius
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma

''''
#A6#
def beta_tin(covalent_radius):##### Really?
    a = 2 * covalent_radius ##### Really ?
    b = a
    c = 2.9 * covalent_radius ##### Really? 
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma
'''

#A-type lattices#
a_lattices = [fcc,bcc,hcp,diamond,bct]




# B-type lattices
#------------------------------------------------------------------------------------------

#B1
def rocksalt(shannon_radius):
    print shannon_radius
    limiting_factors=[2*2**0.2*shannon_radius[0],2*2**0.2*shannon_radius[1],2*shannon_radius[0]+2*shannon_radius[1]]
    a = max(limiting_factors)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma

#B2    
def b2(shannon_radius):
    print shannon_radius
    limiting_factors=[2*(shannon_radius[0]+shannon_radius[0])/np.sqrt(3),2*shannon_radius[1],2*shannon_radius[0]]
    a = max(limiting_factors)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma
    
#B3    
def zincblende(shannon_radius):
    print shannon_radius
    limiting_factors=[2*(max(shannon_radius)*np.sqrt(2)), 4*(shannon_radius[0] + shannon_radius[1])**(1./3.)]
    a = max(limiting_factors)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma


## Zn-S-Zn angle is ~109.5 degrees (from a tetrahedron). It is exactly 2*invCos(-1/3).
## The distance of that Zn-Zn (diagonally to half the face) is (using the cosine rule) is
# root[2(r1+r2)^2 - 2(r1+r2)^(2)cos(ZnSZn angle)].

#B4    
def wurtzite(shannon_radius):
    print shannon_radius
    ##limiting_factors=[2*(max(shannon_radius)*np.sqrt(2)), 4*(shannon_radius[0] + shannon_radius[1])**(1./3.)]
    a = 2*0.827*(shannon_radius[0]+shannon_radius[1])  # 0.817 is sin(109.6/2)
    b = a
    c = (shannon_radius[0]+shannon_radius[1])*(2+2*0.334)  # 0.334 is sin(109.5-90)
    alpha = 90
    beta = 90
    gamma = 120
    return a,b,c,alpha,beta,gamma

#B10
def b10(shannon_radius): #Litharge
    print shannon_radius
    limiting_factors=[4*(max(shannon_radius))/np.sqrt(2), sum(shannon_radius)*1.31]## Explained below.
    a = max(limiting_factors)
    b = a
    c = a*1.26 #Value taken for PbO http://www.mindat.org/min-2466.html#
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,alpha,beta,gamma
## Angle (w) between face-centered atom measured to be ~117.8 deg for   SnO and PbO, where atoms are not similar in size.
b_lattices = [rocksalt,b2,zincblende,wurtzite,b10]

# E-type lattices
#------------------------------------------------------------------------------------------
#E2_1

def E2_1(shannon_radius): #Cubic Pervoskite
    print shannon_radius
    limiting_factors=[2*sum(shannon_radius[1:])]
    a = max(limiting_factors)
    b = a
    c = a
    space = a * np.sqrt(3) - 2 * shannon_radius[1]
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,space,alpha,beta,gamma
e_lattices = [E2_1]


#------------------------------------------------------------------------------------------

def main():

    crystal_elements = raw_input('Which elements? Separate chemical symbols with a space. ')
    crystal_elements = crystal_elements.split()

    perov=False
    if len(crystal_elements) == 3:
        perov = query_yes_no("Are you looking for Perovskites?")
    print perov



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
        print (' \n===============\n ') + lattice.__name__.upper() + (' Lattice \n===============\n  ')
        for elements in crystal_elements:
            x = smact.Element(elements)
            a,b,c,alpha,beta,gamma = lattice(x.covalent_radius)
            print elements,"%.2f"%a,"%.2f"%b,"%.2f"%c,"%.1f"%alpha,"%.1f"%beta,"%.1f"%gamma
    print ' \n Key :: Element a b c alpha beta gamma \n///////////////////////////////////////////////////////////////////////////////////////\n '







    import csv
    oxidation =[]
    coordination =[]

    # Raw input processing

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
        print elements
        print poss_ox_clean
        if len(poss_ox_clean) == 1:
            oxidation.append(int(poss_ox_clean[0]))
        else:
            oxidation.append(int(raw_input('Oxidation state of ' + elements + ' ')))
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
        print elements
        print poss_co_clean
        if len(poss_co_clean) == 1:
            coordination.append(poss_co_clean[0])
        else: 
            coordination.append(raw_input('Coordination state of ' + elements + ' '))

        print elements,oxidation[i],coordination[i]
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
        print (' \n===============\n ') + lattice.__name__.upper() + (' Lattice \n===============\n  ')
        shannon_radius = []
        for i, elements in enumerate(crystal_elements):
            print elements,oxidation[i],coordination[i]
            x = smact.Species(elements,oxidation[i],coordination[i])
            shannon_radius.append(float(x.shannon_radius))
        a,b,c,alpha,beta,gamma = lattice(shannon_radius)
        if lattice.__name__ == 'wurtzite':
            inner_space = a * (6**0.5) - (4*shannon_radius[0])
            print ("With a gap in the middle of " + "%.2f"%inner_space+" (diameter)")
        print crystal_elements[0:2],"%.2f"%a,"%.2f"%b,"%.2f"%c,"%.1f"%alpha,"%.1f"%beta,"%.1f"%gamma



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
        print (' \n===============\n ') + lattice.__name__.upper() + (' Lattice \n===============\n  ')
        shannon_radius = []
        for i, elements in enumerate(crystal_elements):
            print elements,oxidation[i],coordination[i]
            x = smact.Species(elements,oxidation[i],coordination[i])
            shannon_radius.append(float(x.shannon_radius))
        a,b,c,space,alpha,beta,gamma = lattice(shannon_radius)
        print crystal_elements[0:],"%.2f"%a,"%.2f"%b,"%.2f"%c,"%.1f"%alpha,"%.1f"%beta,"%.1f"%gamma
        print ("With a gap in the middle of " + "%.2f"%space+" (diameter)")

if __name__ == '__main__':
    main()
