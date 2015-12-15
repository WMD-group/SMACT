#!/usr/bin/env python

################################################################################
#  Copyright Adam J. Jackson, Daniel Davies (2013)                             #
#                                                                              #
#  This file is part of SMACT: core.py is free software: you can               #
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
import os # get correct path for datafiles when called from another directory
smact_directory = os.path.dirname(__file__)
import csv
# Append a trailing slash to make coherent directory name - this would select the
# root directory in the case of no prefix, so we need to check
if smact_directory:
    smact_directory = smact_directory + '/'

class Element(object):
    """Class providing standard chemical data for elements."""
    def __init__(self, symbol):
        """
        Collection of standard elemental properties for given element.
        Data is drawn from "data/element.txt", part of the Open Babel package.
        element = element('Symbol')

        Atoms with a defined oxidation state draw properties from the "Species" class.

        Attributes:
                
            Element.symbol: Elemental symbol used to retrieve data

            Element.name: Full name of element

            Element.number: Proton number of element

            Element.pauling_eneg: Pauling electronegativity (0.0 if unknown)

            Element.ionpot: Ionisation potential in eV (0.0 if unknown)

            Element.e_affinity: Electron affinity in eV (0.0 if unknown)

            Element.eig: Electron eigenvalue (units unknown)
            N.B. For Cu, Au and Ag this defaults to d-orbital.
            TODO(DD): implement a way to select the desired orbital
            
            Element.crustal_abundance: crustal abundance in the earths crust mg/kg taken from CRC
            
            Element.coord_envs: The allowed coordination enviroments for the ion.

        Raises:
            NameError: Element not found in element.txt
            Warning: Element not found in Eigenvalues.csv
            
        """
        # Import oxidation states from oxidation_states.txt
        with open(smact_directory + 'data/oxidation_states.txt','r') as f:
            data = f.readlines()
        oxidation_states = False
        for line in data:
                        if not line.startswith('#'):
                                l = line.split()
                                if (l[0] == symbol):
                                        if len(l)>1:
                                                oxidation_states = l[1:]
                                                for location, value in enumerate(oxidation_states):
                                                        oxidation_states[location] = int(value)
                                        else:
                                            oxidation_states = []
                                        
        self.oxidation_states = oxidation_states

        # Import crustal abundance from crustal_abundance.txt
        with open(smact_directory + 'data/crustal_abundance.txt','r') as f:
            data = f.readlines()

        crustal_abundance = False

        for line in data:
		if not line.startswith("#"):
		    l = line.split()
		    if (l[0] == symbol):
		 	if len (l)>1:
			    crustal_abundance = l[1]

									
        self.crustal_abundance = crustal_abundance	

        # Import HHI's from HHI's.txt
        with open(smact_directory + "data/HHIs.txt",'r') as f:
	        data = f.readlines() 
        HHI_p, HHI_R = [False, False]
        for line in data:
			if not line.startswith("#"):
				l = line.split()
				if (l[0] == symbol):
					if len(l)>1:
					    HHI_p, HHI_R = l[1:3]
        self.HHI_p = HHI_p
        self.HHI_R = HHI_R								
                        

                                       

        # Import general data from Openbabel-derived data table:
        # Import whole file
        with open(smact_directory + 'data/element.txt','r') as f:
            data = f.readlines()
        # Iterate through data file, ignoring comments and checking line against symbol
        for line in data:
            if not line.startswith('#'):
                l = line.split()
                if (l[1] == symbol):
                    elementdata = l
                    break
        # If the for loop exits without breaking, the element was not found. Report error:
        #else:  
            #raise NameError('Element {0} not found in element.txt'.format(symbol))
        
        # Set attributes
        self.symbol=          str(symbol)
        self.name=            str(elementdata[14])
        self.number=          int(elementdata[0])
        self.covalent_radius= float(elementdata[3])
#        self.vdw_radius=      float(elementdata[5])
        self.pauling_eneg=    float(elementdata[8])
        self.ionpot=          float(elementdata[9])
        self.e_affinity =     float(elementdata[10])

        # Load eigenvalue data from data table by iterating through CSV file
        with open(smact_directory + 'data/Eigenvalues.csv','r') as f:            
            while True:
                l=f.readline()
                if l.split(",")[0] == symbol:
                    self.eig = float(l.split(",")[1])
                    break
                # Check for end of file
                elif not l:
                    #print 'WARNING: Element {0} not found in Eigenvalues.csv'.format(symbol)
                    self.eig = False
                    break

        # Load s-eigenvalue data from data table by iterating through CSV file
        with open(smact_directory + 'data/Eigenvalues_s.csv','r') as f:
            while True:
                l=f.readline()
                if l.split(",")[0] == symbol:
                    self.eig_s = float(l.split(",")[1])
                    break
                # Check for end of file
                elif not l:
                    #print 'WARNING: Element {0} not found in Eigenvalues_s.csv'.format(symbol)
                    self.eig_s = False
                    break

    
#------------------------------------------------------------------------------------------#
# THIS WAS CREATED BY TIM GAUNTLETT. PLEASE CHECK                                          #
#------------------------------------------------------------------------------------------#
        # Load coordination environments
        coord_envs = []
        with open(smact_directory + 'data/shannon_radii.csv','rU') as f:
            reader = csv.reader(f) #, delimiter=',', quoting=csv.QUOTE_NONE)
            for row in reader:
                    if row[0] == symbol:
                        coord_envs.append(row[2])

        self.coord_envs = coord_envs

# Read in the Solid State Energies from csv file
	with open(smact_directory + 'data/SSE.csv','rU') as f: 
            reader = csv.reader(f)#, delimiter=',', quoting=csv.QUOTE_NONE)
            for row in reader:
	        if row[0] == symbol:
	            self.SSE = float(row[2])
	            break
	        elif not l:
	            #print 'WARNING: Element {0} not found in SSE.csv'.format(symbol)
	            self.SSE = False
	            break

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

        # # Load ionic radii data from data table by iterating through CSV file
        # with open(smact_directory + 'data/ionic_radii.csv','r') as f:
        #     while True:
        #         l=f.readline()
        #         if l.split(",")[0] == symbol:
        #             self.ionic = float(l.split(",")[1])
        #             break
        #         # Check for end of file
        #         elif not l:
        #             print 'WARNING: Element {0} not found in ionic_radii.csv'.format(symbol)
        #             self.ionic = False
        #             break

#------------------------------------------------------------------------------------------

class Species(Element):
    """
    Class providing data for elements in a given chemical environment

    In addition to the standard properties from the periodic table (inherited from the 
    Element class), Species objects use the oxidation state and coordination environment
    to provide further properties.

    Attributes: 
        Species.symbol: Elemental symbol used to retrieve data

        Species.name: Full name of element

        Species.oxidation: Oxidation state of species (signed integer)

        Species.coordination: Coordination number of species (integer)

        Species.pauling_eneg: Pauling electronegativity (0.0 if unknown)

        Species.ionpot: Ionisation potential in eV (0.0 if unknown)

        Species.e_affinity: Electron affinity in eV (0.0 if unknown)

        Species.eig: Electron eigenvalue (units unknown)
        N.B. For Cu, Au and Ag this defaults to d-orbital.
        TODO(DD): implement a way to select the desired orbital

    Raises:
        NameError: Element not found in element.txt
        Warning: Element not found in Eigenvalues.csv              

    """

    def __init__(self,symbol,oxidation,coordination):
        Element.__init__(self,symbol)
        self.oxidation = oxidation
        self.coordination = coordination
        # Load shannon radii
        with open(smact_directory + 'data/shannon_radii.csv','rU') as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == symbol and int(row[1]) == oxidation and row[2] == coordination:
                    shannon_radius = row[3]

        self.shannon_radius = shannon_radius


#------------------------------------------------------------------------------------------
def ordered_elements(x,y):
    """ 
Function to return a list of element symbols, ordered by proton number in the range x -> y
Args:
	x,y : integers
Returns:
	ordered_elements : List of element symbols
    """
    with open(smact_directory + 'data/ordered_periodic.txt','r') as f:
    	data = f.readlines()
    elements = []
    for line in data:
		inp = line.split()
		elements.append(inp[0])
	
    ordered_elements = []	
    for i in range(x,y):
        ordered_elements.append(elements[i-1])
		
    return ordered_elements
		


#------------------------------------------------------------------------------------------

def are_eq(A,B,tolerance=1e-4):
    """Check two arrays for tolerance [1,2,3]==[1,2,3]; but [1,3,2]!=[1,2,3]
        Args:
        A/B: arrays
        tolerance: numerical precision for equality condition
        Returns:
        True/False
    """
    are_eq = True
    if len(A) != len(B):
        are_eq = False
    else:
        i = 0
        while i < len(A):
            if abs(A[i] - B[i]) > tolerance:
                are_eq = False
            i = i + 1
    return are_eq

#------------------------------------------------------------------------------------------
def lattices_are_same(lattice1, lattice2):
    """Checks for the equivalence of two lattices
        Args:
        lattice1,lattice2 : ASE crystal class
        Returns:
        lattices_are_same : boolean 
        """
    lattices_are_same = False
    i = 0
    for site1 in lattice1:
        for site2 in lattice2:
            if site1.symbol == site2.symbol:
                if are_eq(site1.position, site2.position):
                    i += 1
    if i == len(lattice1):
        lattices_are_same = True
    return lattices_are_same

#------------------------------------------------------------------------------------------
def charge_neutrality(oxidations, stoichs=False, threshold = 5):
    '''Given a list of oxidation states of arbitrary length it serches for neutral ratios in a given ratio of sites (stoichs) or up to a given threshold
       For the moment it just iterates through binaries, ternaries ... TODO: make it general (KTB)
    Args:
	oxidations : list of integers
	stoichs : stoichiometric ratios for each site (if provided)
	threshold : single threshold to go up to if stoichs are not provided
    Returns:
	exists : bool to say if a ratio exists
	allowed_ratios : ratio that gives neutrality
    '''
    allowed_ratios = []
    ratio_exists = False
    if not stoichs:
	stoichs = []
	for i in range(len(oxidations)):
	    stoichs.append(range(1,threshold + 1))

    if len(oxidations) == 2:
	for i in stoichs[0]:
	    for j in stoichs[1]:	
		if i*oxidations[0] + j*oxidations[1] == 0:
		    if j == 1:
			allowed_ratios.append([i,j])
			ratio_exists = True
# Test to make sure that ratios do not reduce to some simpler ratio
		    elif i%j != 0 and j%1 != 0:
			allowed_ratios.append([i,j])
			ratio_exists = True

    if len(oxidations) == 3:
	for i in stoichs[0]:
	    for j in stoichs[1]:	
	        for k in stoichs[2]:	
		    if i*oxidations[0] + j*oxidations[1] + k*oxidations[2] == 0:
		        if j == 1 and k == 1:
			    allowed_ratios.append([i,j,k])
			    ratio_exists = True
			else:
			    for x in i, j, k:
                                    for y in i, j, k:
                                        if x%y != 0:
			                    allowed_ratios.append([i,j,k])
			        	    ratio_exists = True

    if len(oxidations) == 4:
	for i in stoichs[0]:
	    for j in stoichs[1]:
	    	for k in stoichs[2]:
	    	    for l in stoichs[3]:
			if i*oxidations[0] + j*oxidations[1] + k*oxidations[2] + l*oxidations[3] == 0:
			    if j == 1 and k == 1 and l == 1:
				allowed_ratios.append([i,j,k,l])
				ratio_exists = True
			    else:
			    	for x in i, j, k, l:
				    for y in i, j, k, l:
					if x%y != 0:
					    allowed_ratios.append([i,j,k,l])
                            		    ratio_exists = True	
    if len(oxidations) == 5:
	for i in stoichs[0]:
	    for j in stoichs[1]:
	    	for k in stoichs[2]:
	    	    for l in stoichs[3]:
	    	        for m in stoichs[4]:
			    if i*oxidations[0] + j*oxidations[1] + k*oxidations[2] + l*oxidations[3] + m*oxidations[4] == 0:
			        if j == 1 and k == 1 and l == 1 and m ==1:
				    allowed_ratios.append([i,j,k,l,m])
				    ratio_exists = True
			        else:
			    	    for x in i, j, k, l, m:
				        for y in i, j, k, l, m:
					    if x%y != 0:
					        allowed_ratios.append([i,j,k,l,m])
                            		        ratio_exists = True	



#    return ratio_exists, allowed_ratios
   
#---DWD Alternative output to avoid duplicate ratios---------------------------------------
    unique_ratios = []
    for i in allowed_ratios:
        if i not in unique_ratios:
            unique_ratios.append(i)
    
    return ratio_exists, unique_ratios

 
#------------------------------------------------------------------------------------------
def pauling_test(ox, paul, threshold=0.5):
	''' Testting if a combination of ions makes chemical sense,
	ie positive ions should be of lower Pauling electronegativity
	Args:
	    ox : a list of the oxidation states of the compound
	    paul : the Pauling electronegativities of the elements in the compound
	    threshold : a tolerance for the allowd deviation from the Pauling criterion
	Returns:
	    makes_sense : bool of whether the combination makes sense
	'''
	#return sorted(zip(paul,ox), key=lambda s: s[1])==sorted(zip(paul,ox), key=lambda s: s[0], reverse=True)
        positive = []
	negative = []
	for i, state in enumerate(ox):
	    if state > 0:
		positive.append(paul[i])
	    if state < 0:
		negative.append(paul[i]) 
	if len(positive) == 0 or len(negative) == 0:
	    return False
	if max(positive) - min(negative) > threshold:
	    return False
	else:
	    return True

 
	
#------------------------------------------------------------------------------------------
