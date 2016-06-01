###############################################################################
# Copyright Adam J. Jackson, Daniel Davies (2013)                             #
#                                                                             #
# This file is part of SMACT: smact.__init__ is free software: you can        #
# redistribute it and/or modify it under the terms of the GNU General Public  #
# License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version.                         #
# This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   #
# more details.                                                               #
# You should have received a copy of the GNU General Public License along with#
# this program.  If not, see <http://www.gnu.org/licenses/>.                  #
###############################################################################
"""
Semiconducting Materials from Analogy and Chemical Theory

A collection of fast screening tools from elemental data
"""

# get correct path for datafiles when called from another directory
from os import path
module_directory = path.abspath(path.dirname(__file__))
data_directory = path.join(module_directory, 'data')
import itertools
from itertools import imap

from fractions import gcd
from operator import mul as multiply

from smact import data_loader

class Element(object):
    """Collection of standard elemental properties for given element.

    Data is drawn from "data/element.txt", part of the Open Babel
    package.

    Atoms with a defined oxidation state draw properties from the
    "Species" class.

    Attributes:
        Element.symbol (string): Elemental symbol used to retrieve data

        Element.name (string): Full name of element

        Element.number (int): Proton number of element

        Element.pauling_eneg (float): Pauling electronegativity (0.0 if unknown)

        Element.ionpot (float): Ionisation potential in eV (0.0 if unknown)

        Element.e_affinity (float): Electron affinity in eV (0.0 if unknown)

        Element.eig (float): Electron eigenvalue (units unknown). N.B. For Cu, Au and Ag this defaults to d-orbital.

        Element.eig_s (float): Eigenvalue of s-orbital

        Element.crustal_abundance (float): crustal abundance in the earths crust mg/kg taken from CRC

        Element.coord_envs (list): The allowed coordination enviroments for the ion.

        Element.mass (float) : Molar mass of the element.

        Element.HHI_p (float) : Herfindahl-Hirschman Index for elemental production

        Element.HHI_R (float) : Hirfindahl-Hirschman Index for elemental reserves

        Element.SSE (float) : Solid State Energy

        Element.SSEPauling (float) : SSE based on regression fit with Pauling electronegativity

    Raises:
        NameError: Element not found in element.txt
        Warning: Element not found in Eigenvalues.csv

    """
    def __init__(self, symbol):
        """Initialise Element class

        Args:
            symbol (str): Chemical element symbol (e.g. 'Fe')

        """

        dataset = data_loader.lookup_element_data(symbol, copy=False)

        if dataset == None:
            raise NameError("Elemental data for {0} not found.".format(symbol))

        # Set coordination-environment data from the Shannon-radius data.
        # As above, it is safe to use copy = False with this Get* function.

        shannon_data = data_loader.lookup_element_shannon_radius_data(symbol, copy=False)

        if shannon_data != None:
            coord_envs = [row['coordination'] for row in shannon_data]
        else:
            coord_envs = None

        HHI_scores = data_loader.lookup_element_hhis(symbol)
        if HHI_scores == None:
            HHI_scores = (None, None)

        sse_data = data_loader.lookup_element_sse_data(symbol)
        if sse_data:
            sse = sse_data['SolidStateEnergy']
        else:
            sse = None

        sse_Pauling_data = data_loader.lookup_element_sse_pauling_data(symbol)
        if sse_Pauling_data:
            sse_Pauling = sse_Pauling_data['SolidStateEnergyPauling']
        else:
            sse_Pauling = None

        for attribute, value in (
            ('coord_envs', coord_envs),
            ('covalent_radius', dataset['r_cov']),
            ('crustal_abundance', dataset['Abundance']),
            ('e_affinity', dataset['e_affinity']),
            ('eig', dataset['p_eig']),
            ('eig_s', dataset['s_eig']),
            ('HHI_p', HHI_scores[0]),
            ('HHI_R', HHI_scores[1]),
            ('ionpot', dataset['ion_pot']),
            ('mass', dataset['Mass']),
            ('name', dataset['Name']),
            ('number', dataset['Z']),
            ('oxidation_states',
             data_loader.lookup_element_oxidation_states(symbol)),
            ('pauling_eneg', dataset['el_neg']),
            ('SSE', sse),
            ('SSEPauling', sse_Pauling),
            ('symbol', symbol),
            #('vdw_radius', dataset['RVdW']),
            ):
            setattr(self, attribute, value)

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

    Raises:
        NameError: Element not found in element.txt
        Warning: Element not found in Eigenvalues.csv

    """

    def __init__(self,symbol,oxidation,coordination):
        Element.__init__(self,symbol)

        self.oxidation = oxidation
        self.coordination = coordination

        # Get shannon radius for the oxidation state and coordination.

        self.shannon_radius = None;

        shannon_data = data_loader.lookup_element_shannon_radius_data(symbol);

        for dataset in shannon_data:
            if dataset['charge'] == oxidation and dataset['coordination'] == coordination:
                self.shannon_radius = dataset['crystal_radius'];

        # Get SSE_2015 (revised) for the oxidation state.

        self.SSE_2015 = None

        sse_2015_data = data_loader.lookup_element_sse2015_data(symbol);

        for dataset in sse_2015_data:
            if dataset['OxidationState'] == oxidation:
                self.SSE_2015 = dataset['SolidStateEnergy2015']


def ordered_elements(x,y):
    """
    Return a list of element symbols, ordered by proton number in the range x -> y
    Args:
        x,y : integers
    Returns:
        list: Ordered list of element symbols
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

def element_dictionary(elements=None):
    """
    Create a dictionary of initialised smact.Element objects

    Accessing an Element from a dict is significantly faster than
    repeadedly initialising them on-demand within nested loops.

    Args:
        elements (iterable of strings) : Elements to include. If None,
            use all elements up to 103.
    Returns:
        dict: Dictionary with element symbols as keys and smact.Element
            objects as data
    """
    if elements == None:
        elements = ordered_elements(1,103)
    return {symbol: Element(symbol) for symbol in elements}


def are_eq(A,B,tolerance=1e-4):
    """Check two arrays for tolerance [1,2,3]==[1,2,3]; but [1,3,2]!=[1,2,3]
    Args:
        A, B (lists): 1-D list of values for approximate equality comparison
        tolerance: numerical precision for equality condition

    Returns:
        boolean
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
        boolean
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

    Args:
        oxidations (tuple): Oxidation states of a set of oxidised elements
        stoichs (tuple): Stoichiometry values corresponding to `oxidations`
    """
    return 0 == sum(imap(multiply, oxidations, stoichs))

def neutral_ratios_iter(oxidations, stoichs=False, threshold=5):
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

    Yields:
        tuple: ratio that gives neutrality
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

def neutral_ratios(oxidations, stoichs=False, threshold=5):
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
        threshold (int): Maximum stoichiometry coefficient; if no 'stoichs'
            argument is provided, all combinations of integer coefficients up
            to this value will be tried.

    Returns:
        (exists, allowed_ratios) (tuple):

        exists *bool*:
            True ifc any ratio exists, otherwise False

        allowed_ratios *list of tuples*:
            Ratios of atoms in given oxidation
            states which yield a charge-neutral structure
    """
    allowed_ratios = [x for x in neutral_ratios_iter(oxidations,
                                                        stoichs=stoichs,
                                                        threshold=threshold)]
    return (len(allowed_ratios) > 0, allowed_ratios)
