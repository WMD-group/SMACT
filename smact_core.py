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
import os # get correct path for datafiles when called from another directory
smact_directory = os.path.dirname(__file__)
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

        Attributes:
                
            Element.symbol: Elemental symbol used to retrieve data

            Element.name: Full name of element

            Element.covalent_radius: Covalent radius in AA (1.6 if unknown)
                
            Element.vdw_radius: van der Waals radius in AA (2.0 if unknown)

            Element.pauling_eneg: Pauling electronegativity (0.0 if unknown)

            Element.ionpot: Ionisation potential in eV (0.0 if unknown)

            Element.e_affinity: Electron affinity in eV (0.0 if unknown)

            Element.eig: Electron eigenvalue (units unknown)
            N.B. For Cu, Au and Ag this defaults to d-orbital.
            TODO(DD): implement a way to select the desired orbital

        Raises:
            NameError: Element not found in element.txt
            Warning: Element not found in Eigenvalues.csv
            
        """
        # Import general data from Openbabel-derived data table:
        # Import whole file
        with open(smact_directory + '/data/element.txt','r') as f:
            data = f.readlines()
        # Iterate through data file, ignoring comments and checking line against symbol
        for line in data:
            if not line.startswith('#'):
                l = line.split()
                if (l[1] == symbol):
                    elementdata = l
                    break
        # If the for loop exits without breaking, the element was not found. Report error:
        else:  
            raise NameError('Element {0} not found in element.txt'.format(symbol))
        
        # Set attributes
        self.symbol=          str(symbol)
        self.name=            str(elementdata[14])
        self.covalent_radius= float(elementdata[3])
        self.vdw_radius=      float(elementdata[5])
        self.pauling_eneg=    float(elementdata[8])
        self.ionpot=          float(elementdata[9])
        self.e_affinity =     float(elementdata[10])

        # Load eigenvalue data from data table by iterating through CSV file
        with open(smact_directory + '/data/Eigenvalues.csv','r') as f:            

            while True:
                l=f.readline()
                if l.split(",")[0] == symbol:
                    self.eig = float(l.split(",")[1])
                    break
                # Check for end of file
                elif not l:
                    print 'WARNING: Element {0} not found in Eigenvalues.csv'.format(symbol)
                    self.eig = False
                    break

	# Load s-eigenvalue data from data table by iterating through CSV file
	with open('data/Eigenvalues_s.csv','r') as f:
            while True:
                l=f.readline()
                if l.split(",")[0] == symbol:
                    self.eig_s = float(l.split(",")[1])
                    break
                # Check for end of file
                elif not l:
                    print 'WARNING: Element {0} not found in Eigenvalues_s.csv'.format(symbol)
                    self.eig_s = False
                    break

	# Load ionic radii data from data table by iterating through CSV file
        with open('data/ionic_radii.csv','r') as f:
            while True:
                l=f.readline()
                if l.split(",")[0] == symbol:
                    self.ionic = float(l.split(",")[1])
                    break
                # Check for end of file
                elif not l:
                    print 'WARNING: Element {0} not found in ionic_radii.csv'.format(symbol)
                    self.ionic = False
                    break

#------------------------------------------------------------------------------------------
def are_eq(A,B,tolerance=1e-4):
    """Check two arrays for tolearnce [1,2,3]==[1,2,3]; but [1,3,2]!=[1,2,3]
        Args:
        A/B: arrays
        tolerance: numberical precision for equality condidtion
        Returns:
        True/Flase
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
