from typing import List, Optional, Union

from numpy import product, sqrt

import smact


def eneg_mulliken(element: Union[smact.Element, str]) -> float:
    """Get Mulliken electronegativity from the IE and EA.

    Arguments:
        symbol (smact.Element or str): Element object or symbol

    Returns:
        mulliken (float): Mulliken electronegativity

    """
    if type(element) == str:
        element = smact.Element(element)
    elif type(element) != smact.Element:
        raise Exception(f"Unexpected type: {type(element)}")

    mulliken = (element.ionpot + element.e_affinity) / 2.0

    return mulliken


def band_gap_Harrison(
    anion: str,
    cation: str,
    verbose: bool = False,
    distance: Optional[Union[float, str]] = None,
) -> float:
    """
    Estimates the band gap from elemental data.

    The band gap is estimated using the principles outlined in
    Harrison's 1980 work "Electronic Structure and the Properties of
    Solids: The Physics of the Chemical Bond".

    Args:
        Anion (str): Element symbol of the dominant anion in the system

        Cation (str): Element symbol of the the dominant cation in the system
        Distance (float or str): Nuclear separation between anion and cation
                i.e. sum of ionic radii
        verbose (bool) : An optional True/False flag. If True, additional
        information is printed to the standard output. [Defult: False]

    Returns :
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
    alpha_m = (1.11 * V1_bar) / sqrt(V2**2 + V3**2)

    # Calculate Band gap [(3-43) Harrison 1980 ]
    Band_gap = (3.60 / 3.0) * (sqrt(V2**2 + V3**2)) * (1 - alpha_m)
    if verbose:
        print("V1_bar = ", V1_bar)
        print("V2 = ", V2)
        print("alpha_m = ", alpha_m)
        print("V3 = ", V3)

    return Band_gap


def compound_electroneg(
    verbose: bool = False,
    elements: List[Union[str, smact.Element]] = None,
    stoichs: List[Union[int, float]] = None,
    source: str = "Mulliken",
) -> float:
    """Estimate electronegativity of compound from elemental data.

    Uses Mulliken electronegativity by default, which uses elemental
    ionisation potentials and electron affinities. Alternatively, can
    use Pauling electronegativity, re-scaled by factor 2.86 to achieve
    same scale as Mulliken method (Nethercot, 1974)
    DOI:10.1103/PhysRevLett.33.1088 .

    Geometric mean is used (n-th root of product of components), e.g.:

    X_Cu2S = (X_Cu * X_Cu * C_S)^(1/3)

    Args:
        elements (list) : Elements given as standard elemental symbols.
        stoichs (list) : Stoichiometries, given as integers or floats.
        verbose (bool) : An optional True/False flag. If True, additional information
            is printed to the standard output. [Default: False]
        source: String 'Mulliken' or 'Pauling'; type of Electronegativity to
            use. Note that in SMACT, Pauling electronegativities are
            rescaled to a Mulliken-like scale.

    Returns:
        Electronegativity (float) : Estimated electronegativity (no units).

    """
    if type(elements[0]) == str:
        elementlist = [smact.Element(i) for i in elements]
    elif type(elements[0]) == smact.Element:
        elementlist = elements
    else:
        raise Exception(
            "Please supply a list of element symbols or SMACT Element objects"
        )

    stoichslist = stoichs
    # Convert stoichslist from string to float
    stoichslist = list(map(float, stoichslist))

    # Get electronegativity values for each element

    if source == "Mulliken":
        elementlist = [(el.ionpot + el.e_affinity) / 2.0 for el in elementlist]

    elif source == "Pauling":
        elementlist = [(2.86 * el.pauling_eneg) for el in elementlist]
    else:
        raise Exception(
            f"Electronegativity type '{source}'", "is not recognised"
        )

    # Print optional list of element electronegativities.
    # This may be a useful sanity check in case of a suspicious result.
    if verbose:
        print("Electronegativities of elements=", elementlist)

    # Raise each electronegativity to its appropriate power
    # to account for stoichiometry.
    for i in range(0, len(elementlist)):
        elementlist[i] = [elementlist[i] ** stoichslist[i]]

    # Calculate geometric mean (n-th root of product)
    prod = product(elementlist)
    compelectroneg = (prod) ** (1.0 / (sum(stoichslist)))

    if verbose:
        print("Geometric mean = Compound 'electronegativity'=", compelectroneg)

    return compelectroneg
