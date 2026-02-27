"""
Semiconducting Materials from Analogy and Chemical Theory.

A collection of fast screening tools from elemental data
"""

from __future__ import annotations

__version__ = "4.0.0"

__all__ = [
    "Element",
    "Species",
    "__version__",
    "anions",
    "are_eq",
    "d_block",
    "element_dictionary",
    "lattices_are_same",
    "metals",
    "neutral_ratios",
    "neutral_ratios_iter",
    "ordered_elements",
]

import itertools
import warnings
from functools import reduce
from math import gcd
from operator import mul as multiply
from os import path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

import pandas as pd

module_directory = path.abspath(path.dirname(__file__))
data_directory = path.join(module_directory, "data")
# get correct path for datafiles when called from another directory
from smact import data_loader  # noqa: E402


class Element:
    """
    Collection of standard elemental properties for given element.

    Data is drawn from "data/element.txt", part of the Open Babel
    package.

    Atoms with a defined oxidation state draw properties from the
    "Species" class.

    Attributes:
    ----------
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

        Element.oxidation_states (list) : Default list of allowed oxidation states for use in SMACT. In >3.0, these are the ICSD24 set. In <3.0, these are the SMACT14 set.

        Element.oxidation_states_smact14 (list): Original list of oxidation states that were manually compiled for SMACT in 2014 (default in SMACT < 3.0)

        Element.oxidation_states_sp (list) : List of oxidation states recognised by the Pymatgen Structure Predictor

        Element.oxidation_states_icsd16 (list) : List of oxidation states that appear in the 2016 version of ICSD

        Element.oxidation_states_wiki (list): List of oxidation states that appear wikipedia (https://en.wikipedia.org/wiki/Template:List_of_oxidation_states_of_the_elements) Data retrieved: 2022-09-22

        Element.oxidation_states_custom (list | None ): List of oxidation states that appear in the custom data file supplied (if any)

        Element.oxidation_states_icsd24 (list): List of oxidation states that appear in the 2024 version of the ICSD

        Element.coord_envs (list): The allowed coordination environments for the ion

        Element.covalent_radius (float) : Covalent radius of the element

        Element.mass (float) : Molar mass of the element

        Element.crustal_abundance (float) : Crustal abundance in the earths crust mg/kg taken from CRC

        Element.HHI_p (float) : Herfindahl-Hirschman Index for elemental production

        Element.HHI_r (float) : Hirfindahl-Hirschman Index for elemental reserves

        Element.mendeleev (int): Mendeleev number

        Element.AtomicWeight (float): Atomic weight

        Element.MeltingT (float): Melting temperature in K

        Element.num_valence (int): Number of valence electrons

        Element.num_valence_modified (int): Number of valence electrons based on a modified definition



    Raises:
    ------
        NameError: Element not found in element.txt
        Warning: Element not found in Eigenvalues.csv

    """

    # Instance attribute annotations for static type checking
    coord_envs: list[str] | None
    covalent_radius: float
    crustal_abundance: float
    e_affinity: float
    eig: float
    eig_s: float
    HHI_p: float | None
    HHI_r: float | None
    ionpot: float
    mass: float
    name: str
    number: int
    oxidation_states: list[int] | None
    oxidation_states_smact14: list[int] | None
    oxidation_states_icsd16: list[int] | None
    oxidation_states_sp: list[int] | None
    oxidation_states_wiki: list[int] | None
    oxidation_states_icsd24: list[int] | None
    oxidation_states_custom: list[int] | None
    dipol: float
    pauling_eneg: float | None
    SSE: float | None
    SSEPauling: float | None
    symbol: str
    mendeleev: int | None
    AtomicWeight: float | None
    MeltingT: float | None
    num_valence: int | None
    num_valence_modified: int | None

    def __init__(self, symbol: str, oxi_states_custom_filepath: str | None = None) -> None:
        """
        Initialise Element class.

        Args:
        ----
            symbol (str): Chemical element symbol (e.g. 'Fe')
            oxi_states_custom_filepath (str): Path to custom oxidation states file

        """
        # Get the oxidation states from the custom file if it exists
        if oxi_states_custom_filepath:
            try:
                self._oxidation_states_custom = data_loader.lookup_element_oxidation_states_custom(
                    symbol, oxi_states_custom_filepath
                )
                self.oxidation_states_custom = self._oxidation_states_custom  # type: ignore[assignment]
            except TypeError:
                warnings.warn("Custom oxidation states file not found. Please check the file path.")
                self.oxidation_states_custom = None
        else:
            self.oxidation_states_custom = None
        self.symbol = symbol

        dataset = data_loader.lookup_element_data(self.symbol, copy=False)

        if dataset is None:
            raise NameError(f"Elemental data for {symbol} not found.")

        # Set coordination-environment data from the Shannon-radius data.
        # As above, it is safe to use copy = False with this Get* function.

        shannon_data = data_loader.lookup_element_shannon_radius_data(symbol, copy=False)

        coord_envs = [row["coordination"] for row in shannon_data] if shannon_data is not None else None

        HHI_scores = data_loader.lookup_element_hhis(symbol)
        if HHI_scores is None:
            HHI_scores = (None, None)

        sse_data = data_loader.lookup_element_sse_data(symbol)
        sse = sse_data["SolidStateEnergy"] if sse_data else None

        sse_Pauling_data = data_loader.lookup_element_sse_pauling_data(symbol)
        sse_Pauling = sse_Pauling_data["SolidStateEnergyPauling"] if sse_Pauling_data else None

        magpie_data = data_loader.lookup_element_magpie_data(symbol)
        if magpie_data:
            mendeleev = magpie_data["MendeleevNumber"]
            AtomicWeight = magpie_data["AtomicWeight"]
            MeltingT = magpie_data["MeltingT"]
            num_valence = magpie_data["NValence"]
        else:
            mendeleev = None
            AtomicWeight = None
            MeltingT = None
            num_valence = None

        valence_data = data_loader.lookup_element_valence_data(symbol)
        num_valence_modified = valence_data["NValence"] if valence_data else None

        self.coord_envs = coord_envs
        self.covalent_radius = dataset["r_cov"]
        self.crustal_abundance = dataset["Abundance"]
        self.e_affinity = dataset["e_affinity"]
        self.eig = dataset["p_eig"]
        self.eig_s = dataset["s_eig"]
        self.HHI_p = HHI_scores[0]
        self.HHI_r = HHI_scores[1]
        self.ionpot = dataset["ion_pot"]
        self.mass = dataset["Mass"]
        self.name = dataset["Name"]
        self.number = dataset["Z"]
        self.oxidation_states = data_loader.lookup_element_oxidation_states_icsd24(symbol)
        self.oxidation_states_smact14 = data_loader.lookup_element_oxidation_states(symbol)
        self.oxidation_states_icsd16 = data_loader.lookup_element_oxidation_states_icsd(symbol)
        self.oxidation_states_sp = data_loader.lookup_element_oxidation_states_sp(symbol)
        self.oxidation_states_wiki = data_loader.lookup_element_oxidation_states_wiki(symbol)
        self.oxidation_states_icsd24 = data_loader.lookup_element_oxidation_states_icsd24(symbol)
        self.dipol = dataset["dipol"]
        self.pauling_eneg = dataset["el_neg"]
        self.SSE = sse
        self.SSEPauling = sse_Pauling
        self.mendeleev = mendeleev
        self.AtomicWeight = AtomicWeight
        self.MeltingT = MeltingT
        self.num_valence = num_valence
        self.num_valence_modified = num_valence_modified


class Species(Element):
    """
    Class providing data for elements in a given chemical environment.

    In addition to the standard properties from the periodic table
    (inherited from the  Element class), Species objects use the
    oxidation state and coordination environment to provide further
    properties.
    The Species object can be created with either a default set of shannon radii (radii_source='shannon') or with a set of machine-learnt shannon radii (radii_source='extended').
    The source of the machine-learnt shannon radii set is
    Baloch, A.A., Alqahtani, S.M., Mumtaz, F., Muqaibel, A.H., Rashkeev, S.N. and Alharbi, F.H., 2021. Extending Shannon's ionic radii database using machine learning. Physical Review Materials, 5(4), p.043804.

    Attributes:
    ----------
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
    ------
        NameError: Element not found in element.txt
        Warning: Element not found in Eigenvalues.csv

    """

    def __init__(
        self,
        symbol: str,
        oxidation: int,
        coordination: int = 4,
        radii_source: str = "shannon",
    ) -> None:
        """
        Initialise Species class.

        Args:
        ----
        symbol (str): Chemical element symbol (e.g. 'Fe')
        oxidation (int): Oxidation state of species
        coordination (int): Coordination number of species
        radii_source (str): Source of shannon radii data. Choose 'shannon' for
            the default shannon radii set or 'extended' for the machine-learnt shannon radii set

        """
        Element.__init__(self, symbol)

        self.oxidation = oxidation
        self.coordination = coordination

        # Get shannon radius for the oxidation state and coordination.

        self.shannon_radius = None
        self.ionic_radius = None

        if radii_source == "shannon":
            shannon_data = data_loader.lookup_element_shannon_radius_data(symbol)
        elif radii_source == "extended":
            shannon_data = data_loader.lookup_element_shannon_radius_data_extendedML(symbol)
        else:
            raise ValueError(f"Data source {radii_source!r} not recognised. Choose 'shannon' or 'extended'.")

        if shannon_data:
            for dataset in shannon_data:
                if dataset["charge"] == oxidation and str(coordination) == dataset["coordination"].split("_")[0]:
                    self.shannon_radius = dataset["crystal_radius"]
                    self.ionic_radius = dataset["ionic_radius"]
                    break

        # Get the average shannon and ionic radii
        self.average_shannon_radius = None
        self.average_ionic_radius = None

        if shannon_data:
            # Get the rows of the shannon radius table for the element
            shannon_data_df = pd.DataFrame(shannon_data)

            # Get the rows corresponding to the oxidation state of the species
            charge_rows = shannon_data_df.loc[shannon_data_df["charge"] == oxidation]

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


def ordered_elements(x: int, y: int) -> list[str]:
    """
    Return a list of element symbols, ordered by proton number in the range x -> y
    Args:
        x,y : integers
    Returns:
        list: Ordered list of element symbols.
    """
    with open(path.join(data_directory, "ordered_periodic.txt")) as f:
        data = f.readlines()
    elements = [line.split()[0] for line in data]

    return elements[x - 1 : y]


def element_dictionary(
    elements: Iterable[str] | None = None,
    oxi_states_custom_filepath: str | None = None,
) -> dict[str, Element]:
    """
    Create a dictionary of initialised smact.Element objects.

    Accessing an Element from a dict is significantly faster than
    repeadedly initialising them on-demand within nested loops.

    Args:
    ----
        elements (iterable of strings) : Elements to include. If None,
            use all elements up to 103.
        oxi_states_custom_filepath (str): Path to custom oxidation states file


    Returns:
    -------
        dict: Dictionary with element symbols as keys and smact.Element
            objects as data

    """
    if elements is None:
        elements = ordered_elements(1, 103)
    if oxi_states_custom_filepath:
        return {symbol: Element(symbol, oxi_states_custom_filepath) for symbol in elements}
    else:
        return {symbol: Element(symbol) for symbol in elements}


def are_eq(A: Sequence[float], B: Sequence[float], tolerance: float = 1e-4) -> bool:
    """Check two arrays for approximate element-wise equality.

    Args:
        A: 1-D sequence of values.
        B: 1-D sequence of values.
        tolerance: Absolute tolerance for equality.

    Returns:
        True if arrays are element-wise equal within tolerance, False otherwise.
    """
    if len(A) != len(B):
        return False
    return bool(np.allclose(A, B, atol=tolerance, rtol=0))


def lattices_are_same(lattice1: Sequence, lattice2: Sequence, tolerance: float = 1e-4) -> bool:
    """Check whether two ASE lattices contain the same sites.

    Args:
        lattice1: ASE crystal class.
        lattice2: ASE crystal class.
        tolerance: Absolute tolerance for position equality.

    Returns:
        True if every site in *lattice1* has a matching site in *lattice2*.
    """
    matches = 0
    for site1 in lattice1:
        for site2 in lattice2:
            if site1.symbol == site2.symbol and are_eq(site1.position, site2.position, tolerance=tolerance):
                matches += 1
    return matches == len(lattice1)


def _gcd_recursive(*args: int) -> int:
    """Get the greatest common denominator among any number of ints."""
    return reduce(gcd, args)


def _isneutral(oxidations: Sequence[int], stoichs: Sequence[int]) -> bool:
    """
    Check if set of oxidation states is neutral in given stoichiometry.

    Args:
    ----
        oxidations (tuple): Oxidation states of a set of oxidised elements
        stoichs (tuple): Stoichiometry values corresponding to `oxidations`

    """
    return sum(map(multiply, oxidations, stoichs)) == 0


def neutral_ratios_iter(
    oxidations: Sequence[int],
    stoichs: Sequence[Sequence[int]] | None = None,
    threshold: int | None = 5,
) -> filter:
    """
    Iterator for charge-neutral stoichiometries.

    Given a list of oxidation states of arbitrary length, yield ratios in which
    these form a charge-neutral compound. Stoichiometries may be provided as a
    set of legal stoichiometries per site (e.g. a known family of compounds);
    otherwise all unique ratios are tried up to a threshold coefficient.

    Args:
    ----
        oxidations : list of integers
        stoichs : stoichiometric ratios for each site (if provided)
        threshold : single threshold to go up to if stoichs are not provided

    Yields:
    ------
        tuple: ratio that gives neutrality

    """
    if stoichs is None:
        if threshold is None:
            raise ValueError("threshold must be an int when stoichs is not provided")
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
    oxidations: Sequence[int],
    stoichs: Sequence[Sequence[int]] | None = None,
    threshold: int | None = 5,
) -> list[tuple[int, ...]]:
    """
    Get a list of charge-neutral compounds.

    Given a list of oxidation states of arbitrary length, return ratios in which
    these form a charge-neutral compound. Stoichiometries may be provided as a
    set of legal stoichiometries per site (e.g. a known family of compounds);
    otherwise all unique ratios are tried up to a threshold coefficient.

    Args:
    ----
        oxidations (list of ints): Oxidation state of each site
        stoichs (list of positive ints): A selection of valid stoichiometric
            ratios for each site
        threshold (int): Maximum stoichiometry coefficient; if no 'stoichs'
            argument is provided, all combinations of integer coefficients up
            to this value will be tried.

    Returns:
    -------
        list of tuples: Ratios of atoms in given oxidation states which yield
            a charge-neutral structure. Empty list if no ratios exist.

    """
    return list(neutral_ratios_iter(oxidations, stoichs=stoichs, threshold=threshold))


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
