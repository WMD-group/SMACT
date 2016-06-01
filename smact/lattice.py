#!/usr/bin/env python 

################################################################################
#  Copyright Keith T. Butler, Adam J. Jackson (2013)                           #
#                                                                              #
# This file is part of SMACT: lattice.py is free software: you can             #
# redistribute it and/or modify it under the terms of the GNU General Public   #
# License as published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.                          #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for    #
# more details.                                                                #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
################################################################################

import numpy as np
import smact

class Lattice(object):
      """A unique set of Sites

      Lattice objects define a general crystal structure, with a space group and
      a collection of Site objects. These Site objects have their own fractional
      coordinates and a list of possible oxidation states (see the Site class).

      Specific crystal structures with elements assigned to sites are
      "materials" and use the Atoms class from the Atomic Simulation
      Environment.

      Attributes: 
          basis_sites: A list of Site objects [SiteA, SiteB, SiteC, ...]
          comprising the basis sites in Cartesian coordinates

          space_group: Integer space group number according to the
          International Tables for Crystallography.  

          structurbericht:
          Structurbericht identity, if applicable (e.g. 'B1')

      Methods:
          lattice_vector_calc():

      """

      def __init__(self, sites, space_group=1, strukturbericht=False):
            self.sites = sites
            self.space_group = space_group
            self.strukturbericht = strukturbericht

class Site(object):
      """
      A single lattice site with a list of possible oxidation states

      The Site object is primarily used within Lattice objects.

      Attributes:
          position: A list of fractional coordinates [x,y,z]
          oxidation_states: A list of possible oxidation states e.g. [-1,0,1]
      
      """ 
      
      def __init__(self, position, oxidation_states=[0]):
            self.position = position
            self.oxidation_states = oxidation_states


########## Everything below this is probably broken ##########

#------------------------------------------------------------------------------------
def check_lattice_charges(charges, site_elements, sites):
      """
      This function checks the sum of the charges on the lattice sites in the crystal.
      The formal oxidation states of each species are assumed.
      
      'This does not signify that the chemical bonds in the crystal 
      are necessarily ionic in the sense of the quantum mechanics'
      Linus Pauling (1929)

      Args:
          charges: array of the charge contributions of each lattice site
          site_elements: All possible elements which satisfy charge neutrality 
          sites: list, the current set of elements under inspection

      Returns:
          input site_elements array with appended neutral sites (???)

      """
# Check if all sub-lattice charges sum to zero, append to list if true
      if np.sum(charges) == 0:
          site_elements.append(sites)
# Need to use two lists, new contains the present composition; site_elements is the 
# list of all news found.

      return site_elements
#------------------------------------------------------------------------------------
def possible_compositions(crystal, elements):
    """
    Search for the elements which satisfy the possible oxidation states and 
    provide charge neutrality.

    Args:
        crystal: A Lattice object defining the crystal class
        elements: Dictionary of elements and their allowable oxidation states

    Returns:
        list/array/dict of int/float/string of something useful (???)

    """

    for site in crystal.sites:
	composition = []
	for ox in site.oxidation_states:
	    total_charge = total_charge + ox


'''
# Initialise the array atom, containing possible elements for each sub lattice
    atom = []
# Initialise the array, site_elements, containing compositions found
    site_elements = []
    i = 0
    while i < len(crystal.sites):
        atom.append(possible_elements(elements, crystal.sites[i].oxidation_states))
        i = i + 1
    print atom
for site in crystal
# I could not think of an elegant way to generalise this. For now it loops through the possible
# lattices until it reaches the number, then it goes no further. Therefore, many of the loops
# here are redundant. 
# site_n refers to the possible element @ site n
# we multiply these by oxidation number to get sub-lattice charge
    i = 0
    for site_1 in atom[0]:
        charges = np.zeros(shape=(len(atom)))
# Look up the charge of 'site_1' multiply by the multiplicity of that site to get charge
        charges[0] = int(elements[site_1]) * crystal.site_ratios[0]
	if len(atom) == 1: #Are ther more sites to check?
	    sites = [site_1] 
	    site_elements = check_lattice_charges(charges, site_elements, sites)
	if 2 <= len(atom):       #Are ther more sites to check?
      	    for site_2 in atom[1]:
                charges[1] = int(elements[site_2]) * crystal.site_ratios[1]
		if len(atom) == 2:    #If we have sampled all sites test the charge
	            sites = [site_1, site_2]
		    site_elements = check_lattice_charges(charges, site_elements, sites)
	        if 3 <= len(atom):    #Are ther more sites to check?
                    for site_3 in atom[2]:
                        charges[2] = int(elements[site_3]) * crystal.site_ratios[2]
	                if len(atom) == 3: #If we have sampled all sites test the charge
	                    sites = [site_1, site_2, site_3]
		            site_elements = check_lattice_charges(charges, site_elements, sites)
		        if 4 <= len(atom): #Are ther more sites to check?
                            for site_4 in atom[3]:
                                charges[3] = int(elements[site_4]) * crystal.site_ratios[3]
	                        if len(atom) == 4: 
	                            sites = [site_1, site_2, site_3, site_4]
		                    site_elements = check_lattice_charges(charges, site_elements, sites)
			    if 5 <= len(atom):
			    	for site_5 in atom[4]:
				    charges[4] = int(elements[site_5]) * crystal.site_ratios[4] 
	                            if len(atom) == 5:
	                                sites = [site_1, site_2, site_3, site_4]
		                        site_elements = check_lattice_charges(charges, site_elements, sites)

    return site_elements
'''
#------------------------------------------------------------------------------------
def possible_elements(elements, oxidations):
    """Identify possible atoms to occupy a site
    
    Args:
        elements:
        oxidations: 

    Returns: 
        Array of atoms

    """
    atoms = []
    for element in elements:
	elemental_oxidations = smact.Element(element).oxidation_states
        for ox_state_a in oxidations:
	    for element_ox in elemental_oxidations:
        	if int(element_ox) == int(ox_state_a):
        	    atoms.append(element)
    return atoms
#------------------------------------------------------------------------------------
