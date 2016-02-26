################################################################################
#  Copyright Adam J. Jackson, Daniel Davies (2013)                             #
#                                                                              #
#  This file is part of SMACT: smact.__init__ is free software: you can        #
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
Semiconducting Materials from Analogy and Chemical Theory

A collection of fast screening tools from elemental data
"""

from os import path # get correct path for datafiles when called from another directory
module_directory = path.abspath(path.dirname(__file__))
data_directory = path.join(module_directory, 'data')
import itertools

import csv
from fractions import gcd
from operator import mul as multiply

class Element(object):
    """Class providing standard chemical data for elements."""
    def __init__(self, symbol):
        """
        Collection of standard elemental properties for given element.
        Data is drawn from "data/element.txt", part of the Open Babel package.
        element = element('Symbol')

        Atoms with a defined oxidation state draw properties from the "Species" class.

        Attributes:
            Element.symbol (string): Elemental symbol used to retrieve data

            Element.name (string): Full name of element

            Element.number (int): Proton number of element

            Element.pauling_eneg (float): Pauling electronegativity (0.0 if unknown)

            Element.ionpot (float): Ionisation potential in eV (0.0 if unknown)

            Element.e_affinity (float): Electron affinity in eV (0.0 if unknown)

            Element.eig (float): Electron eigenvalue (units unknown)
            N.B. For Cu, Au and Ag this defaults to d-orbital.
            TODO(DD): implement a way to select the desired orbital

            Element.eig_s (float): Eigenvalue of s-orbital

            Element.crustal_abundance (float): crustal abundance in the earths crust mg/kg taken from CRC

            Element.coord_envs (list): The allowed coordination enviroments for the ion.

        Raises:
            NameError: Element not found in element.txt
            Warning: Element not found in Eigenvalues.csv

        """
        # Import oxidation states from oxidation_states.txt
        with open(path.join(data_directory,
                            'oxidation_states.txt'),'r') as f:
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
        with open(path.join(data_directory,
                            'crustal_abundance.txt'), 'r') as f:
            data = f.readlines()

        crustal_abundance = False

        for line in data:
            if not line.startswith("#"):
                l = line.split()
                if (l[0] == symbol) and len(l) > 1:
                    crustal_abundance = float(l[1])

        self.crustal_abundance = crustal_abundance

        # Import HHI's from HHI's.txt
        with open(path.join(data_directory,
                            "HHIs.txt"), 'r') as f:
            data = f.readlines()
        HHI_p, HHI_R = [False, False]
        for line in data:
            if not line.startswith("#"):
                l = line.split()
                if (l[0] == symbol):
                    if len(l)>1:
                        HHI_p, HHI_R = l[1:3]
        self.HHI_p = float(HHI_p)
        self.HHI_R = float(HHI_R)

        # Import general data from Openbabel-derived data table:
        # Import whole file
        with open(path.join(data_directory,
                            'element.txt'), 'r') as f:
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
        with open(path.join(data_directory,
                            'Eigenvalues.csv'),'r') as f:
            while True:
                l=f.readline()
                if l.split(",")[0] == symbol:
                    self.eig = float(l.split(",")[1])
                    break
                # Check for end of file
                elif not l:
                    self.eig = False
                    break

        # Load s-eigenvalue data from data table by iterating through CSV file
        with open(path.join(data_directory,
                            'Eigenvalues_s.csv'), 'r') as f:
            while True:
                l=f.readline()
                if l.split(",")[0] == symbol:
                    self.eig_s = float(l.split(",")[1])
                    break
                # Check for end of file
                elif not l:
                    self.eig_s = False
                    break


#----------------------------------------------------------------------------#
# THIS WAS CREATED BY TIM GAUNTLETT. PLEASE CHECK                            #
#----------------------------------------------------------------------------#
        # Load coordination environments
        coord_envs = []
        with open(path.join(data_directory,
                            'shannon_radii.csv'), 'rU') as f:
            reader = csv.reader(f)
            for row in reader:
                    if row[0] == symbol:
                        coord_envs.append(row[2])

        self.coord_envs = coord_envs

# Read in the Solid State Energies from csv file
        with open(path.join(data_directory, 'SSE.csv'), 'rU') as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == symbol:
                    self.SSE = float(row[2])
                    break
                elif not l:
                    self.SSE = False
                    break

#------------------------------------------------------------------------------------------

class Species(Element):
    """
    Class providing data for elements in a given chemical environment

    In addition to the standard properties from the periodic table
    (inherited from the  Element class), Species objects use the
    oxidation state and coordination environment to provide further
    properties.

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
        with open(path.join(data_directory,
                            'shannon_radii.csv'),'rU') as f:
            reader = csv.reader(f)
            for row in reader:
                if (row[0] == symbol
                    and int(row[1]) == oxidation
                    and row[2] == coordination):
                    shannon_radius = row[3]

        self.shannon_radius = shannon_radius


def ordered_elements(x,y):
    """
    Return a list of element symbols, ordered by proton number in the range x -> y
    Args:
        x,y : integers
    Returns:
        ordered_elements : List of element symbols
    """
    with open(path.join(data_directory,
                        'ordered_periodic.txt'), 'r') as f:
        data = f.readlines()
    elements = []
    for line in data:
        inp = line.split()
        elements.append(inp[0])

    ordered_elements = []
    for i in range(x,y+1):
        ordered_elements.append(elements[i-1])

    return ordered_elements

def element_dictionary(elements):
    '''
    Given a list of element names this returns a dictionary containing all of the information about those elements. This can lead to significant speedups in iterative scripts where the element object is constantly referred to.
    Args:
        elements : a list of the element names
    Returns:
        element_dictionary : a dictionary with element names as keys and smact.Elements as data,
    '''
    dictionary = {}
    for ele in elements:
        dictionary[ele] = Element(ele)
    return dictionary
    

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


def lattices_are_same(lattice1, lattice2, tolerance=1e-4):
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
                if are_eq(site1.position,
                          site2.position,
                          tolerance=tolerance):
                    i += 1
    if i == len(lattice1):
        lattices_are_same = True
    return lattices_are_same


def _gcd_recursive(*args):
    """
    Get the greatest common denominator among any number of ints
    """
    if len(args) == 2:
        return gcd(*args)
    else:
        return gcd(args[0], _gcd_recursive(*args[1:]))

def _isneutral(oxidations, stoichs):
    """
    Check if set of oxidation states is neutral in given stoichiometry
    """
    return 0 == sum(itertools.imap(multiply, oxidations, stoichs))
    
def charge_neutrality_iter(oxidations, stoichs=False, threshold=5):
    """
    Iterator for charge-neutral stoichiometries
    
    Given a list of oxidation states of arbitrary length, yield ratios in which
    these form a charge-neutral compound. Stoichiometries may be provided as a
    set of legal stoichiometries per site (e.g. a known family of compounds);
    otherwise all unique ratios are tried up to a threshold coefficient.

    Args:
        oxidations : list of integers
        stoichs : stoichiometric ratios for each site (if provided)
        threshold : single threshold to go up to if stoichs are not provided

    Returns:
        exists : bool to say if a ratio exists
        allowed_ratios : ratio that gives neutrality
    """    
    if not stoichs:
        stoichs = [range(1,threshold+1)] * len(oxidations)

    # First filter: remove combinations which have a common denominator
    # greater than 1 (i.e. Use simplest form of each set of ratios)
    # Second filter: return only charge-neutral combinations    
    return itertools.ifilter(
        lambda x: _isneutral(oxidations, x) and _gcd_recursive(*x) == 1,
        # Generator: enumerate all combinations of stoichiometry
        itertools.product(*stoichs)
        )

def charge_neutrality(oxidations, stoichs=False, threshold=5):
    """
    Get a list of charge-neutral compounds
    
    Given a list of oxidation states of arbitrary length, yield ratios in which
    these form a charge-neutral compound. Stoichiometries may be provided as a
    set of legal stoichiometries per site (e.g. a known family of compounds);
    otherwise all unique ratios are tried up to a threshold coefficient.

    Given a list of oxidation states of arbitrary length it searches for
    neutral ratios in a given ratio of sites (stoichs) or up to a given
    threshold.

    Args:
        oxidations (list of ints): Oxidation state of each site
        stoichs (list of positive ints): A selection of valid stoichiometric
            ratios for each site
        threshold (int): Maximum stoichiometry coefficient; if no "stoichs"
            argument is provided, all combinations of integer coefficients up
            to this value will be tried.

    Returns:
        exists (bool): True if any ratio exists, otherwise False
        allowed_ratios (list of tuples): Ratios of atoms in given oxidation
            states which yield a charge-neutral structure
    """
    allowed_ratios = [x for x in charge_neutrality_iter(oxidations,
                                                        stoichs=stoichs,
                                                        threshold=threshold)]
    return (len(allowed_ratios) > 0,
            allowed_ratios)


def pauling_test(ox, paul, threshold=0.5):
    """ Testing if a combination of ions makes chemical sense,
    ie positive ions should be of lower Pauling electronegativity

    Args:
        ox : a list of the oxidation states of the compound
        paul : the Pauling electronegativities of the elements in the compound
        threshold : a tolerance for the allowed deviation from the Pauling
        criterion

    Returns:
        makes_sense : bool of whether the combination makes sense

    """
    positive = []
    negative = []
    for i, state in enumerate(ox):
        if state > 0:
            positive.append(paul[i])
        if state < 0:
	    if paul[i] in negative: # Reject materials where the same anion occupies two sites.
                return False
            else:
                negative.append(paul[i])
    if len(positive) == 0 or len(negative) == 0:
        return False
    if max(positive) == min(negative):
        return False
    if max(positive) - min(negative) > threshold:
        return False
    else:
        return True
