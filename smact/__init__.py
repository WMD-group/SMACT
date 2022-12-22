"""
Semiconducting Materials from Analogy and Chemical Theory

A collection of fast screening tools from elemental data
"""

import itertools
import warnings
from math import gcd
from operator import mul as multiply
from os import path
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import pandas as pd

module_directory = path.abspath(path.dirname(__file__))
data_directory = path.join(module_directory, "data")
# get correct path for datafiles when called from another directory
from smact import data_loader


class Element:
    """Collection of standard elemental properties for given element.

    Data is drawn from "data/element.txt", part of the Open Babel
    package.

    Atoms with a defined oxidation state draw properties from the
    "Species" class.

    Attributes:
        Element.symbol (string) : Elemental symbol used to retrieve data

        Element.name (string) : Full name of element

        Element.number (int) : Proton number of element

        Element.pauling_eneg (float) : Pauling electronegativity (0.0 if unknown)

        Element.ionpot (float) : Ionisation potential in eV (0.0 if unknown)

        Element.e_affinity (float) : Electron affinity in eV (0.0 if unknown)

        Element.dipol (float) : Static dipole polarizability in 1.6488e-41 C m^2 / V  (0.0 if unknown)

        Element.eig (float) : Electron eigenvalue (units unknown) N.B. For Cu, Au and Ag this defaults to d-orbital

        Element.eig_s (float) : Eigenvalue of s-orbital

        Element.SSE (float) : Solid State Energy

        Element.SSEPauling (float) : SSE based on regression fit with Pauling electronegativity

        Element.oxidation_states (list) : Default list of allowed oxidation states for use in SMACT

        Element.oxidation_states_sp (list) : List of oxdation states recognised by the Pymatgen Structure Predictor

        Element.oxidation_states_icsd (list) : List of oxidation states that appear in the ICSD

        Element.oxidation_states_wiki (list): List of oxidation states that appear wikipedia (https://en.wikipedia.org/wiki/Template:List_of_oxidation_states_of_the_elements) Data retrieved: 2022-09-22

        Element.coord_envs (list): The allowed coordination enviroments for the ion

        Element.covalent_radius (float) : Covalent radius of the element

        Element.mass (float) : Molar mass of the element

        Element.crustal_abundance (float) : Crustal abundance in the earths crust mg/kg taken from CRC

        Element.HHI_p (float) : Herfindahl-Hirschman Index for elemental production

        Element.HHI_r (float) : Hirfindahl-Hirschman Index for elemental reserves

    Raises:
        NameError: Element not found in element.txt
        Warning: Element not found in Eigenvalues.csv

    """

    def __init__(self, symbol: str):
        """Initialise Element class

        Args:
            symbol (str): Chemical element symbol (e.g. 'Fe')

        """

        dataset = data_loader.lookup_element_data(symbol, copy=False)

        if dataset == None:
            raise NameError(f"Elemental data for {symbol} not found.")

        # Set coordination-environment data from the Shannon-radius data.
        # As above, it is safe to use copy = False with this Get* function.

        shannon_data = data_loader.lookup_element_shannon_radius_data(
            symbol, copy=False
        )

        if shannon_data != None:
            coord_envs = [row["coordination"] for row in shannon_data]
        else:
            coord_envs = None

        HHI_scores = data_loader.lookup_element_hhis(symbol)
        if HHI_scores == None:
            HHI_scores = (None, None)

        sse_data = data_loader.lookup_element_sse_data(symbol)
        if sse_data:
            sse = sse_data["SolidStateEnergy"]
        else:
            sse = None

        sse_Pauling_data = data_loader.lookup_element_sse_pauling_data(symbol)
        if sse_Pauling_data:
            sse_Pauling = sse_Pauling_data["SolidStateEnergyPauling"]
        else:
            sse_Pauling = None

        for attribute, value in (
            ("coord_envs", coord_envs),
            ("covalent_radius", dataset["r_cov"]),
            ("crustal_abundance", dataset["Abundance"]),
            ("e_affinity", dataset["e_affinity"]),
            ("eig", dataset["p_eig"]),
            ("eig_s", dataset["s_eig"]),
            ("HHI_p", HHI_scores[0]),
            ("HHI_r", HHI_scores[1]),
            ("ionpot", dataset["ion_pot"]),
            ("mass", dataset["Mass"]),
            ("name", dataset["Name"]),
            ("number", dataset["Z"]),
            (
                "oxidation_states",
                data_loader.lookup_element_oxidation_states(symbol),
            ),
            (
                "oxidation_states_icsd",
                data_loader.lookup_element_oxidation_states_icsd(symbol),
            ),
            (
                "oxidation_states_sp",
                data_loader.lookup_element_oxidation_states_sp(symbol),
            ),
            (
                "oxidation_states_wiki",
                data_loader.lookup_element_oxidation_states_wiki(symbol),
            ),
            ("dipol", dataset["dipol"]),
            ("pauling_eneg", dataset["el_neg"]),
            ("SSE", sse),
            ("SSEPauling", sse_Pauling),
            ("symbol", symbol),
            # ('vdw_radius', dataset['RVdW']),
        ):
            setattr(self, attribute, value)


class Species(Element):
    """
    Class providing data for elements in a given chemical environment

    In addition to the standard properties from the periodic table
    (inherited from the  Element class), Species objects use the
    oxidation state and coordination environment to provide further
    properties.
    The Species object can be created with either a default set of shannon radii (radii_source='shannon') or with a set of machine-learnt shannon radii (radii_source='extended').
    The source of the machine-learnt shannon radii set is
    Baloch, A.A., Alqahtani, S.M., Mumtaz, F., Muqaibel, A.H., Rashkeev, S.N. and Alharbi, F.H., 2021. Extending Shannon's ionic radii database using machine learning. Physical Review Materials, 5(4), p.043804.

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

        Species.shannon_radius: Shannon radius of Species.

        Species.ionic_radius: Ionic radius of Species.

        Species.average_shannon_radius: An average shannon radius for the species. The average is taken over all coordination environments.

        Species.average_ionic_radius: An average ionic radius for the species. The average is taken over all coordination environments.

    Raises:
        NameError: Element not found in element.txt
        Warning: Element not found in Eigenvalues.csv

    """

    def __init__(
        self,
        symbol: str,
        oxidation: int,
        coordination: int = 4,
        radii_source: str = "shannon",
    ):
        Element.__init__(self, symbol)

        self.oxidation = oxidation
        self.coordination = coordination

        # Get shannon radius for the oxidation state and coordination.

        self.shannon_radius = None

        if radii_source == "shannon":

            shannon_data = data_loader.lookup_element_shannon_radius_data(
                symbol
            )

        elif radii_source == "extended":
            shannon_data = (
                data_loader.lookup_element_shannon_radius_data_extendedML(
                    symbol
                )
            )

        else:
            print(
                "Data source not recognised. Please select 'shannon' or 'extended'. "
            )

        if shannon_data:
            for dataset in shannon_data:
                if (
                    dataset["charge"] == oxidation
                    and str(coordination)
                    == dataset["coordination"].split("_")[0]
                ):
                    self.shannon_radius = dataset["crystal_radius"]

        # Get ionic radius
        self.ionic_radius = None

        if shannon_data:
            for dataset in shannon_data:
                if (
                    dataset["charge"] == oxidation
                    and str(coordination)
                    == dataset["coordination"].split("_")[0]
                ):
                    self.ionic_radius = dataset["ionic_radius"]

        # Get the average shannon and ionic radii
        self.average_shannon_radius = None
        self.average_ionic_radius = None

        if shannon_data:
            # Get the rows of the shannon radius table for the element
            shannon_data_df = pd.DataFrame(shannon_data)

            # Get the rows corresponding to the oxidation state of the species
            charge_rows = shannon_data_df.loc[
                shannon_data_df["charge"] == oxidation
            ]

            # Get the mean
            self.average_shannon_radius = charge_rows["crystal_radius"].mean()
            self.average_ionic_radius = charge_rows["ionic_radius"].mean()

        # Get SSE_2015 (revised) for the oxidation state.

        self.SSE_2015 = None

        sse_2015_data = data_loader.lookup_element_sse2015_data(symbol)
        if sse_2015_data:
            for dataset in sse_2015_data:
                if dataset["OxidationState"] == oxidation:
                    self.SSE_2015 = dataset["SolidStateEnergy2015"]
        else:
            self.SSE_2015 = None


def ordered_elements(x: int, y: int) -> List[str]:
    """
    Return a list of element symbols, ordered by proton number in the range x -> y
    Args:
        x,y : integers
    Returns:
        list: Ordered list of element symbols
    """
    with open(path.join(data_directory, "ordered_periodic.txt")) as f:
        data = f.readlines()
    elements = []
    for line in data:
        inp = line.split()
        elements.append(inp[0])

    ordered_elements = []
    for i in range(x, y + 1):
        ordered_elements.append(elements[i - 1])

    return ordered_elements


def element_dictionary(elements: Optional[Iterable[str]] = None):
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
        elements = ordered_elements(1, 103)
    return {symbol: Element(symbol) for symbol in elements}


def are_eq(A: list, B: list, tolerance: float = 1e-4):
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
                if are_eq(site1.position, site2.position, tolerance=tolerance):
                    i += 1
    if i == len(lattice1):
        lattices_are_same = True
    return lattices_are_same


def _gcd_recursive(*args: Iterable[int]):
    """
    Get the greatest common denominator among any number of ints
    """
    if len(args) == 2:
        return gcd(*args)
    else:
        return gcd(args[0], _gcd_recursive(*args[1:]))


def _isneutral(oxidations: Tuple[int, ...], stoichs: Tuple[int, ...]):
    """
    Check if set of oxidation states is neutral in given stoichiometry

    Args:
        oxidations (tuple): Oxidation states of a set of oxidised elements
        stoichs (tuple): Stoichiometry values corresponding to `oxidations`
    """
    return 0 == sum(map(multiply, oxidations, stoichs))


def neutral_ratios_iter(
    oxidations: List[int],
    stoichs: Union[bool, List[int]] = False,
    threshold: Optional[int] = 5,
):
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
        stoichs = [list(range(1, threshold + 1))] * len(oxidations)

    # First filter: remove combinations which have a common denominator
    # greater than 1 (i.e. Use simplest form of each set of ratios)
    # Second filter: return only charge-neutral combinations
    return filter(
        lambda x: _isneutral(oxidations, x) and _gcd_recursive(*x) == 1,
        # Generator: enumerate all combinations of stoichiometry
        itertools.product(*stoichs),
    )


def neutral_ratios(
    oxidations: List[int], stoichs: Union[bool, List[int]] = False, threshold=5
):
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
    allowed_ratios = [
        x
        for x in neutral_ratios_iter(
            oxidations, stoichs=stoichs, threshold=threshold
        )
    ]
    return (len(allowed_ratios) > 0, allowed_ratios)


# List of metals
metals = [
    "Li",
    "Be",
    "Na",
    "Mg",
    "Al",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
]

# List of elements that can be considered 'anions'.
# Similar to the Pymatgen 'electronegative elements' but excluding H, B, C & Si.
anions = ["N", "P", "As", "Sb", "O", "S", "Se", "Te", "F", "Cl", "Br", "I"]

# List of d-block metals
d_block = [
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "La",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
]
