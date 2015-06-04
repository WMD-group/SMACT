#!/usr/bin/env python 
# Using the ase spacegroup module this can build the structure, from 
# the composition, as defined in the smact_lattice module. 
#TO DO:
# Calculating the lattice parameters from the elemnts involved
################################################################################
# Copyright Tim Gauntlett, Adam J. Jackson, Keith Butler   (2014)              #
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

import smact.core as core
import smact.lattice_parameters as lattice_parameters
import numpy as np
import os
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm



# get correct path for datafiles when called from another directory
smact_directory = '../../smact'
# Append a trailing slash to make coherent directory name - this would select the
# root directory in the case of no prefix, so we need to check
if smact_directory:
    smact_directory = smact_directory + '/'


def radius_covalent(element):
    with open(smact_directory + 'data/Covalent_radii.csv','rU') as f:
        reader = csv.reader(f)         
        reader.next() # Skip the first row which is titles   
        r_covalent=0 
        for row in reader:
            symbols= list(row[0])
            # Some symbols specify a shell e.g. "Cu(d)"; we can ignore this part
            for char in symbols:
                if char == '(':
                    a=symbols.index("(")
                    symbols = symbols[:a]
            symbols="".join(symbols)
            if symbols == element:
                r_covalent = float(row[3])/ 100 
    return r_covalent
        
def formal_charge(element):
    """Given an elemental symbol, return a list of available oxidation states"""
    with open(smact_directory + 'data/oxidation_states.txt','r') as f:
        data = f.readlines()
        oxidation_states = False
        for line in data:
            if not line.startswith('#'):
                l = line.split()
                if (l[0] == element):
                    if len(l)>1:
                        oxidation_states = l[1:]
                        for location, value in enumerate(oxidation_states):
                            oxidation_states[location] = int(value)
                    else:
                        oxidation_states = []
    return oxidation_states

def radius_shannon(element,oxidation):
    with open(smact_directory + 'data/shannon_radii.csv','rU') as f:
        reader = csv.reader(f)
        r_shannon=False
        for row in reader:
            if element == row[0] and float(oxidation)==float(row[1]) and row[2]=='4_n': #only for Wurtzite
                r_shannon=row[3]
    return r_shannon                 
    
            
#------------------------------------------------------------------------------------------
output=[]

# Processing Lattices
#------------------------------------------------------------------------------------------
lattice_charge=-1
upper_n = 100
lower_n = 1

all_elements=[]
element_set=core.ordered_elements(lower_n,upper_n)     ## The safety valve!

# Form a separate complete list of elements for quickly accessing proton numbers with index()
with open(smact_directory + 'data/ordered_periodic.txt','rU') as f:
    reader = csv.reader(f)
    for row in reader:
        all_elements.append(row[0])

# Form list of elements known to coordinate tetragonally
tetra_coord_elements=[]
for element in element_set:
    with open(smact_directory + 'data/shannon_radii.csv','rU') as f:
        reader = csv.reader(f)
        r_shannon=False
        for row in reader:
            if row[2]=="4_n" and row[0]==element and row[0] not in tetra_coord_elements:
                tetra_coord_elements.append(row[0])

# Calculate properties of combinations with target charge
holes=np.zeros(shape=(100,100))
matrix_hhi=np.zeros(shape=(100,100))

for a_element in tetra_coord_elements:
    a_cov=radius_covalent(a_element)
    a_formal_charge = formal_charge(a_element)
    b_elements = [b for b in tetra_coord_elements if (
                        tetra_coord_elements.index(b) > tetra_coord_elements.index(a_element))]
    if a_formal_charge:#i.e. a_formal_charge has a value
        for a_ox in a_formal_charge:
            a_shan=radius_shannon(a_element,a_ox)

            if a_shan == False:  # Fall back to covalent radius if Shannon unavailable
                for b_element in b_elements:
                    b_formal_charge = formal_charge(b_element)  ## check again
                    if b_formal_charge: #i.e. b_formal_charge has a value
                        # Check for charge-neutral combinations and write to array
                        # Note that only the elements are recorded; 
                        # multiple charge combinations for one element pair will be overwritten
                        for b_ox in (x for x in b_formal_charge if x+a_ox==lattice_charge):
                            b_cov=radius_covalent(b_element)
                            g=[float(a_cov),float(b_cov)]
                            method="covalent radii"
                            a,b,c,alpha,beta,gamma = lattice_parameters.wurtzite(g)
			# Calculate the space inside the cage
			    inner_space = inner_space = a * (6**0.5) - (4*g[0])
                            hhi_a=core.Element(a_element).HHI_p
                            hhi_b=core.Element(b_element).HHI_p
                            hhi_rp=(float(hhi_a)*float(hhi_b))**0.5
                            holes[all_elements.index(a_element),all_elements.index(b_element)]=inner_space
                            matrix_hhi[all_elements.index(a_element),all_elements.index(b_element)]=hhi_rp
                            output.append([a_element,b_element,hhi_rp,method,a,b,c,inner_space])
                            
            else:    # Use Shannon radius
                for b_element in b_elements:
                    b_formal_charge = formal_charge(b_element)  ## check again
                    if b_formal_charge: 
                        # Check for charge-neutral combinations and write to array
                        for b_ox in b_formal_charge:
                            if radius_shannon(b_element,b_ox) == False:
                                if a_ox+b_ox==lattice_charge:                       ## Charge of lattice
                                    b_cov=radius_covalent(b_element)
                                    g=[float(a_cov),float(b_cov)]
                                    method="covalent radii"
                                    a,b,c,alpha,beta,gamma = lattice_parameters.wurtzite(g)
			    # Calculate the space inside the cage
			            inner_space = inner_space = a * (6**0.5) - (4*g[0])
                                    hhi_a=core.Element(a_element).HHI_p
                                    hhi_b=core.Element(b_element).HHI_p
                                    hhi_rp=(float(hhi_a)*float(hhi_b))**0.5
                                    holes[all_elements.index(a_element),all_elements.index(b_element)]=inner_space
                                    matrix_hhi[all_elements.index(a_element),all_elements.index(b_element)]=hhi_rp
                                    output.append([a_element,b_element,hhi_rp,method,a,b,c,inner_space])
                                    
                                    
                            else:   ### Shannon radius
                                b_shan=radius_shannon(b_element,b_ox)
                                if a_ox+b_ox==lattice_charge:
                                    g=[float(a_shan),float(b_shan)]
                                    method="shannon radii"
                                    a,b,c,alpha,beta,gamma = lattice_parameters.wurtzite(g)
			    # Calculate the space inside the cage
			            inner_space = inner_space = a * (6**0.5) - (4*g[0])
                                    hhi_a=core.Element(a_element).HHI_p
                                    hhi_b=core.Element(b_element).HHI_p
                                    hhi_rp=(float(hhi_a)*float(hhi_b))**0.5
                                    holes[all_elements.index(a_element),all_elements.index(b_element)]=inner_space
                                    matrix_hhi[all_elements.index(a_element),all_elements.index(b_element)]=hhi_rp
                                    output.append([a_element,b_element,hhi_rp,method,a,b,c,inner_space])
                                
with open('wurtzite.csv', 'w') as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=',',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writerow(["Element A", "Element B", "HHI index (geometric mean)",
                         "method","a / AA","b / AA","c / AA","pore volume / AA^3"])
    for structure in output:
        csv_writer.writerow(structure)


hhi_max = 2000
hhi_min = 1     # To exclude unlisted elements (e.g. Tc...)

hhi_filtered_output = [compound for compound in output if ((compound[2] <= hhi_max)
                                                           and (compound[2] >= hhi_min))]

def get_volume_key(compound):
    return compound[7]
size_sorted_ouput = sorted(hhi_filtered_output, key=get_volume_key, reverse=True)

with open('candidates.csv','w') as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=',',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writerow(["Element A", "Element B", "HHI index (geometric mean)",
                         "method","a / AA","b / AA","c / AA","pore volume / AA^3"])
    csv_writer.writerows(size_sorted_ouput)


mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['font.size'] = 16

plt.figure()
plt.imshow(matrix_hhi, cmap=cm.BuPu, interpolation='none',
           extent=[lower_n, upper_n, lower_n, upper_n],
           origin='lower')
plt.xlabel('Elements(A) by proton numbers',fontsize=22)
plt.ylabel('Elements(B) by proton numbers',fontsize=22)
plt.xlim(0, upper_n)
plt.ylim(0, upper_n)
cbar = plt.colorbar()
cbar.set_label("Herfindahl-Hirschman Index (geometric mean)",
               rotation=270, labelpad=20)
plt.savefig('wurtzite_hhi.png')
#plt.show()
plt.figure()
'''
plt.imshow(holes, cmap=cm.BuGn, interpolation='none',
           extent=[lower_n, upper_n, lower_n, upper_n],
           origin='lower')
plt.xlabel('Elements(A) by proton numbers',fontsize=22)
plt.ylabel('Elements(B) by proton numbers',fontsize=22)
plt.xlim(0, upper_n)
plt.ylim(0, upper_n)
cbar = plt.colorbar()
cbar.set_label("Est. pore volume / AA^3",
               rotation=270, labelpad=20)
plt.savefig('wurtzite_holes.png')
#plt.show()
'''
