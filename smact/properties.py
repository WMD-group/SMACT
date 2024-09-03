"""A collection of tools for estimating physical properties based on chemical composition."""

from __future__ import annotations

import numpy as np

import smact
from smact.utils.composition import parse_formula


def eneg_mulliken(element: smact.Element | str) -> float:
    """
    Get Mulliken electronegativity from the IE and EA.

    Arguments:
    ---------
        element (smact.Element or str): Element object or symbol

    Returns:
    -------
        mulliken (float): Mulliken electronegativity

    """
    if isinstance(element, str):
        element = smact.Element(element)
    elif not isinstance(element, smact.Element):
        raise TypeError(f"Unexpected type: {type(element)}")

    return (element.ionpot + element.e_affinity) / 2.0


def band_gap_Harrison(
    anion: str,
    cation: str,
    verbose: bool = False,
    distance: float | str | None = None,
) -> float:
    """
    Estimates the band gap from elemental data.

    The band gap is estimated using the principles outlined in
    Harrison's 1980 work "Electronic Structure and the Properties of
    Solids: The Physics of the Chemical Bond".

    Args:
    ----
        anion (str): Element symbol of the dominant anion in the system
        cation (str): Element symbol of the the dominant cation in the system
        distance (float or str): Nuclear separation between anion and cation
                i.e. sum of ionic radii
        verbose (bool) : An optional True/False flag. If True, additional
            information is printed to the standard output. [Default: False]

    Returns:
    -------
        Band_gap (float): Band gap in eV

    """
    # Set constants
    hbarsq_over_m = 7.62

    # Get anion and cation
    An = anion
    Cat = cation
    d = float(distance)

    # Get elemental data:
    elements_dict = smact.element_dictionary((An, Cat))
    An, Cat = elements_dict[An], elements_dict[Cat]

    # Calculate values of equation components
    V1_Cat = (Cat.eig - Cat.eig_s) / 4
    V1_An = (An.eig - An.eig_s) / 4
    V1_bar = (V1_An + V1_Cat) / 2
    V2 = 2.16 * hbarsq_over_m / (d**2)
    V3 = (Cat.eig - An.eig) / 2
    alpha_m = (1.11 * V1_bar) / np.sqrt(V2**2 + V3**2)

    # Calculate Band gap [(3-43) Harrison 1980 ]
    Band_gap = (3.60 / 3.0) * (np.sqrt(V2**2 + V3**2)) * (1 - alpha_m)
    if verbose:
        print("V1_bar = ", V1_bar)
        print("V2 = ", V2)
        print("alpha_m = ", alpha_m)
        print("V3 = ", V3)

    return Band_gap


def compound_electroneg(
    verbose: bool = False,
    elements: list[str | smact.Element] | None = None,
    stoichs: list[int | float] | None = None,
    source: str = "Mulliken",
) -> float:
    """
    Estimate electronegativity of compound from elemental data.

    Uses Mulliken electronegativity by default, which uses elemental
    ionisation potentials and electron affinities. Alternatively, can
    use Pauling electronegativity, re-scaled by factor 2.86 to achieve
    same scale as Mulliken method (Nethercot, 1974)
    DOI:10.1103/PhysRevLett.33.1088 .

    Geometric mean is used (n-th root of product of components), e.g.:

    X_Cu2S = (X_Cu * X_Cu * C_S)^(1/3)

    Args:
    ----
        elements (list) : Elements given as standard elemental symbols.
        stoichs (list) : Stoichiometries, given as integers or floats.
        verbose (bool) : An optional True/False flag. If True, additional information
            is printed to the standard output. [Default: False]
        source: String 'Mulliken' or 'Pauling'; type of Electronegativity to
            use. Note that in SMACT, Pauling electronegativities are
            rescaled to a Mulliken-like scale.

    Returns:
    -------
        Electronegativity (float) : Estimated electronegativity (no units).

    """
    if isinstance(elements[0], str):
        elementlist = [smact.Element(i) for i in elements]
    elif isinstance(elements[0], smact.Element):
        elementlist = elements
    else:
        raise TypeError("Please supply a list of element symbols or SMACT Element objects")

    stoichslist = stoichs
    # Convert stoichslist from string to float
    stoichslist = list(map(float, stoichslist))

    # Get electronegativity values for each element

    if source == "Mulliken":
        elementlist = [(el.ionpot + el.e_affinity) / 2.0 for el in elementlist]

    elif source == "Pauling":
        elementlist = [(2.86 * el.pauling_eneg) for el in elementlist]
    else:
        raise Exception(f"Electronegativity type '{source}'", "is not recognised")

    # Print optional list of element electronegativities.
    # This may be a useful sanity check in case of a suspicious result.
    if verbose:
        print("Electronegativities of elements=", elementlist)

    # Raise each electronegativity to its appropriate power
    # to account for stoichiometry.
    for i in range(len(elementlist)):
        elementlist[i] = [elementlist[i] ** stoichslist[i]]

    # Calculate geometric mean (n-th root of product)
    prod = np.prod(elementlist)
    compelectroneg = (prod) ** (1.0 / (sum(stoichslist)))

    if verbose:
        print("Geometric mean = Compound 'electronegativity'=", compelectroneg)

    return compelectroneg


def valence_electron_count(compound: str) -> float:
    """
    Calculate the Valence Electron Count (VEC) for a given chemical compound.

    This function parses the input compound, extracts the elements and their
    stoichiometries, and calculates the VEC using the valence electron data
    from SMACT's Element class.

    Args:
        compound (str): Chemical formula of the compound (e.g., "Fe2O3").

    Returns:
        float: Valence Electron Count (VEC) for the compound.

    Raises:
        ValueError: If an element in the compound is not found in the valence data.
    """

    def get_element_valence(element: str) -> int:
        try:
            return smact.Element(element).num_valence_modified
        except NameError:
            raise ValueError(f"Valence data not found for element: {element}") from None

    element_stoich = parse_formula(compound)

    total_valence = 0
    total_stoich = 0
    for element, stoich in element_stoich.items():
        try:
            valence = get_element_valence(element)
            total_valence += stoich * valence
            total_stoich += stoich
        except TypeError:
            raise ValueError(f"No valence information for element {element}")

    if total_stoich == 0:
        return 0.0

    return total_valence / total_stoich
