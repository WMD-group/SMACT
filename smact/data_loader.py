"""
Provide data from text files while transparently caching for efficiency.

This module handles the loading of external data used to initialise the
core smact.Element and smact.Species classes.  It implements a
transparent data-caching system to avoid a large amount of I/O when
naively constructing several of these objects.  It also implements a
switchable system to print verbose warning messages about possible
missing data (mainly for debugging purposes). In general these functions
are used in the background and it is not necessary to use them directly.
"""

from __future__ import annotations

import csv
import functools
import os
import warnings

import pandas as pd

from smact import data_directory

__all__ = [
    "float_or_None",
    "lookup_element_data",
    "lookup_element_hhis",
    "lookup_element_magpie_data",
    "lookup_element_oxidation_states",
    "lookup_element_oxidation_states_custom",
    "lookup_element_oxidation_states_icsd",
    "lookup_element_oxidation_states_icsd24",
    "lookup_element_oxidation_states_sp",
    "lookup_element_oxidation_states_wiki",
    "lookup_element_shannon_radius_data",
    "lookup_element_shannon_radius_data_extendedML",
    "lookup_element_sse2015_data",
    "lookup_element_sse_data",
    "lookup_element_sse_pauling_data",
    "lookup_element_valence_data",
    "set_warnings",
]

# Module-level switch: emit verbose warning messages about missing data.
_print_warnings = False


def set_warnings(enable=True):
    """
    Set verbose warning messages on and off.

    In order to see any of the warnings, this function needs to be
    called _before_ the first call to the smact.Element()
    constructor.

    Args:
    ----
    enable (bool) : print verbose warning messages.

    """
    global _print_warnings
    _print_warnings = enable


def _warn(message: str) -> None:
    """Emit a warning through the warnings module when verbose mode is on."""
    if _print_warnings:
        warnings.warn(message, stacklevel=3)


def _get_data_rows(filename):
    """Generator for datafile entries by row for custom oxidation states lists. Skips rows with no oxidation states for performance."""
    with open(filename) as file:
        for line in file:
            line = line.strip()
            if line and line[0] != "#" and any(char.isdigit() for char in line):
                yield line.split()


def float_or_None(x):
    """Cast a string to a float or to a None."""
    try:
        return float(x)
    except ValueError:
        return None


@functools.cache
def _load_oxidation_states(filepath: str) -> dict[str, list[int]]:
    """Load an oxidation-state text file into a symbol→states dict, cached by filepath."""
    data: dict[str, list[int]] = {}
    for items in _get_data_rows(filepath):
        data[items[0]] = [int(x) for x in items[1:]]
    return data


def _lookup_ox_states(
    symbol: str,
    filepath: str,
    copy: bool = True,
) -> list[int] | None:
    """Generic helper: look up oxidation states for *symbol* from *filepath*."""
    data = _load_oxidation_states(filepath)
    if symbol in data:
        return list(data[symbol]) if copy else data[symbol]
    _warn(f"Oxidation states for element {symbol} not found.")
    return None


def lookup_element_oxidation_states(symbol: str, copy: bool = True) -> list[int] | None:
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used is the SMACT default (smact14) and
    most exhaustive list.

    Args:
        symbol (str): the atomic symbol of the element to look up.
        copy (bool): if True (default), return a copy of the oxidation-state
            list rather than a reference to the cached data.

    Returns:
        list: Known oxidation states, or None if not found.
    """
    return _lookup_ox_states(symbol, os.path.join(data_directory, "oxidation_states.txt"), copy)


def lookup_element_oxidation_states_icsd(symbol: str, copy: bool = True) -> list[int] | None:
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those found
    in the ICSD (and judged to be non-spurious).

    Args:
        symbol (str): the atomic symbol of the element to look up.
        copy (bool): if True (default), return a copy of the list.

    Returns:
        list: Known oxidation states, or None if not found.
    """
    return _lookup_ox_states(symbol, os.path.join(data_directory, "oxidation_states_icsd.txt"), copy)


def lookup_element_oxidation_states_sp(symbol: str, copy: bool = True) -> list[int] | None:
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those that
    are in the Pymatgen default lambda table for structure prediction.

    Args:
        symbol (str): the atomic symbol of the element to look up.
        copy (bool): if True (default), return a copy of the list.

    Returns:
        list: Known oxidation states, or None if not found.
    """
    return _lookup_ox_states(symbol, os.path.join(data_directory, "oxidation_states_SP.txt"), copy)


def lookup_element_oxidation_states_wiki(symbol: str, copy: bool = True) -> list[int] | None:
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those that
    appear on Wikipedia (https://en.wikipedia.org/wiki/Template:List_of_oxidation_states_of_the_elements).

    Args:
        symbol (str): the atomic symbol of the element to look up.
        copy (bool): if True (default), return a copy of the list.

    Returns:
        list: Known oxidation states, or None if not found.
    """
    return _lookup_ox_states(symbol, os.path.join(data_directory, "oxidation_states_wiki.txt"), copy)


def lookup_element_oxidation_states_custom(
    symbol: str,
    filepath: str,
    copy: bool = True,
) -> list[int] | list[str] | dict[str, list[int]] | None:
    """
    Retrieve a list of known oxidation states for an element from a user-supplied file.

    The cache is keyed by *filepath*, so calling with different files returns
    the correct data for each file (unlike the old single-global behaviour).

    Args:
        symbol (str): the atomic symbol to look up. Pass ``"all"`` to return
            all symbols in the file.
        filepath (str): path to the text file containing oxidation-state data.
        copy (bool): if True (default), return a copy of the list.

    Returns:
        list: Known oxidation states for the element, or None if not found.
    """
    data = _load_oxidation_states(filepath)
    if symbol == "all":
        return list(data) if copy else data
    if symbol in data:
        return list(data[symbol]) if copy else data[symbol]
    _warn(f"Oxidation states for element {symbol} not found.")
    return None


def lookup_element_oxidation_states_icsd24(symbol: str, copy: bool = True) -> list[int] | None:
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those found
    in the 2024 version of the ICSD (with ≥5 reports).

    Args:
        symbol (str): the atomic symbol of the element to look up.
        copy (bool): if True (default), return a copy of the list.

    Returns:
        list: Known oxidation states, or None if not found.
    """
    return _lookup_ox_states(
        symbol,
        os.path.join(data_directory, "oxidation_states_icsd24_filtered.txt"),
        copy,
    )


@functools.cache
def _load_hhis() -> dict[str, tuple[float, float]]:
    data: dict[str, tuple[float, float]] = {}
    with open(os.path.join(data_directory, "hhi.txt")) as file:
        for line in file:
            line = line.strip()
            if line and line[0] != "#":
                items = line.split()
                data[items[0]] = (float(items[1]), float(items[2]))
    return data


def lookup_element_hhis(symbol):
    """
    Retrieve the HHI_R and HHI_p scores for an element.

    Args:
    ----
        symbol : the atomic symbol of the element to look up.

    Returns:
    -------
        tuple : (HHI_p, HHI_R)

            Return None if values for the elements were
            not found in the external data.

    """
    data = _load_hhis()
    if symbol in data:
        return data[symbol]
    _warn(f"HHI data for element {symbol} not found.")
    return None


_ELEMENT_DATA_KEYS = (
    "Symbol",
    "Name",
    "Z",
    "Mass",
    "r_cov",
    "e_affinity",
    "p_eig",
    "s_eig",
    "Abundance",
    "el_neg",
    "ion_pot",
    "dipol",
)


@functools.cache
def _load_element_data() -> dict[str, dict]:
    data: dict[str, dict] = {}
    for items in _get_data_rows(os.path.join(data_directory, "element_data.txt")):
        clean_items = items[0:2] + list(map(float_or_None, items[2:]))
        data[items[0]] = dict(list(zip(_ELEMENT_DATA_KEYS, clean_items, strict=False)))
    return data


def lookup_element_data(symbol: str, copy: bool = True):
    """
    Retrieve tabulated data for an element.

    The table "data/element_data.txt" contains a collection of relevant
    atomic data.

    Args:
    ----
        symbol (str) : Atomic symbol for lookup
        copy (bool) : if True (default), return a copy of the
            data dictionary, rather than a reference to the cached
            object -- only used copy=False in performance-sensitive code
            and where you are certain the dictionary will not be
            modified!

    Returns:
    -------
        dict: Dictionary of data for given element, keyed by column headings from data/element_data.txt.

    """
    data = _load_element_data()
    if symbol in data:
        return data[symbol].copy() if copy else data[symbol]
    _warn(f"Elemental data for {symbol} not found.")
    return None


@functools.cache
def _load_shannon_radii() -> dict[str, list[dict]]:
    data: dict[str, list[dict]] = {}
    with open(os.path.join(data_directory, "shannon_radii.csv")) as file:
        reader = csv.reader(file)
        next(reader)  # skip header
        for row in reader:
            key = row[0]
            dataset = {
                "charge": int(row[1]),
                "coordination": row[2],
                "crystal_radius": float(row[3]),
                "ionic_radius": float(row[4]),
                "comment": row[5],
            }
            if key in data:
                data[key].append(dataset)
            else:
                data[key] = [dataset]
    return data


def lookup_element_shannon_radius_data(symbol, copy=True):
    """
    Retrieve Shannon radii for known states of an element.

    Retrieve Shannon radii for known oxidation states and coordination
    environments of an element.

    Args:
    ----
        symbol (str) : the atomic symbol of the element to look up.

        copy (Optional(bool)): if True (default), return a copy of the data
        dictionary, rather than a reference to the cached object --
        only use copy=False in performance-sensitive code and where
        you are certain the dictionary will not be modified!

    Returns:
    -------
        list:
            Shannon radii datasets.

        Returns None if the element was not found among the external
        data.

        Shannon radii datasets are dictionaries with the keys:

        charge
            *int* charge
        coordination
            *int* coordination
        crystal_radius
            *float*
        ionic_radius
            *float*
        comment
            *str*

    """
    data = _load_shannon_radii()
    if symbol in data:
        return [item.copy() for item in data[symbol]] if copy else data[symbol]
    _warn(f"Shannon-radius data for element {symbol} not found.")
    return None


@functools.cache
def _load_shannon_radii_extendedML() -> dict[str, list[dict]]:
    data: dict[str, list[dict]] = {}
    with open(os.path.join(data_directory, "shannon_radii_ML_extended.csv")) as file:
        reader = csv.reader(file)
        next(reader)  # skip header
        for row in reader:
            key = row[0]
            dataset = {
                "charge": int(row[1]),
                "coordination": row[2],
                "crystal_radius": float(row[3]),
                "ionic_radius": float(row[4]),
                "comment": row[5],
            }
            if key in data:
                data[key].append(dataset)
            else:
                data[key] = [dataset]
    return data


def lookup_element_shannon_radius_data_extendedML(symbol, copy=True):
    """
    Retrieve the machine learned extended Shannon radii for
    known states of an element.


    Retrieve Shannon radii for known oxidation states and coordination
    environments of an element.

    Source of extended radii is:
    Baloch, A.A., Alqahtani, S.M., Mumtaz, F., Muqaibel, A.H., Rashkeev,
    S.N. and Alharbi, F.H., 2021.
    Extending Shannon's Ionic Radii Database Using Machine Learning.
    arXiv preprint arXiv:2101.00269.

    Args:
    ----
        symbol (str) : the atomic symbol of the element to look up.

        copy (Optional(bool)): if True (default), return a copy of the data
        dictionary, rather than a reference to the cached object --
        only use copy=False in performance-sensitive code and where
        you are certain the dictionary will not be modified!

    Returns:
    -------
        list:
            Extended Shannon radii datasets.

        Returns None if the element was not found among the external
        data.

        Shannon radii datasets are dictionaries with the keys:

        charge
            *int* charge
        coordination
            *int* coordination
        crystal_radius
            *float*
        ionic_radius
            *float*
        comment
            *str*

    """
    data = _load_shannon_radii_extendedML()
    if symbol in data:
        return [item.copy() for item in data[symbol]] if copy else data[symbol]
    _warn(f"Extended Shannon-radius data for element {symbol} not found.")
    return None


@functools.cache
def _load_sse_data() -> dict[str, dict]:
    data: dict[str, dict] = {}
    with open(os.path.join(data_directory, "SSE.csv")) as file:
        reader = csv.reader(file)
        for row in reader:
            data[row[0]] = {
                "AtomicNumber": int(row[1]),
                "SolidStateEnergy": float(row[2]),
                "IonisationPotential": float(row[3]),
                "ElectronAffinity": float(row[4]),
                "MullikenElectronegativity": float(row[5]),
                "SolidStateRenormalisationEnergy": float(row[6]),
            }
    return data


def lookup_element_sse_data(symbol):
    """
    Retrieve the solid-state energy (SSE) data for an element.

    Taken from J. Am. Chem. Soc., 2011, 133 (42), pp 16852-16960,
    DOI: 10.1021/ja204670s

    Args:
    ----
        symbol : the atomic symbol of the element to look up.

    Returns:
    -------
        list : SSE datasets for the element, or None
            if the element was not found among the external data.

        SSE datasets are dictionaries with the keys:

        AtomicNumber
            *int*
        SolidStateEnergy
            *float* SSE
        IonisationPotential
            *float*
        ElectronAffinity
            *float*
        MullikenElectronegativity
            *str*
        SolidStateRenormalisationEnergy
            *float*

    """
    data = _load_sse_data()
    if symbol in data:
        return data[symbol].copy()
    _warn(f"Solid-state energy data for element {symbol} not found.")
    return None


@functools.cache
def _load_sse2015_data() -> dict[str, list[dict]]:
    data: dict[str, list[dict]] = {}
    with open(os.path.join(data_directory, "SSE_2015.csv")) as file:
        reader = csv.reader(file)
        for row in reader:
            key = row[0]
            dataset = {
                "OxidationState": int(row[1]),
                "SolidStateEnergy2015": float(row[2]),
            }
            if key in data:
                data[key].append(dataset)
            else:
                data[key] = [dataset]
    return data


def lookup_element_sse2015_data(symbol, copy=True):
    """
    Retrieve SSE (2015) data for element in oxidation state.

    Retrieve the solid-state energy (SSE2015) data for an element in an
    oxidation state.  Taken from J. Solid State Chem., 2015, 231,
    pp138-144, DOI: 10.1016/j.jssc.2015.07.037.

    Args:
    ----
        symbol : the atomic symbol of the element to look up.
        copy: if True (default), return a copy of the data dictionary,
        rather than a reference to a cached object -- only use
        copy=False in performance-sensitive code and where you are
        certain the dictionary will not be modified!

    Returns:
    -------
        list : SSE datasets for the element, or None
            if the element was not found among the external data.

        SSE datasets are dictionaries with the keys:

        OxidationState
            *int*
        SolidStateEnergy2015
            *float* SSE2015

    """
    data = _load_sse2015_data()
    if symbol in data:
        return [item.copy() for item in data[symbol]] if copy else data[symbol]
    _warn(f"Solid-state energy (revised 2015) data for element {symbol} not found.")
    return None


@functools.cache
def _load_sse_pauling_data() -> dict[str, dict]:
    data: dict[str, dict] = {}
    with open(os.path.join(data_directory, "SSE_Pauling.csv")) as file:
        reader = csv.reader(file)
        for row in reader:
            data[row[0]] = {"SolidStateEnergyPauling": float(row[1])}
    return data


def lookup_element_sse_pauling_data(symbol):
    """
    Retrieve Pauling SSE data.

    Retrieve the solid-state energy (SSEPauling) data for an element
    from the regression fit when SSE2015 is plotted against Pauling
    electronegativity.  Taken from J. Solid State Chem., 2015, 231,
    pp138-144, DOI: 10.1016/j.jssc.2015.07.037

    Args:
    ----
    symbol (str) : the atomic symbol of the element to look up.

    Returns: A dictionary containing the SSE2015 dataset for the
        element, or None if the element was not found among the external
        data.

    """
    data = _load_sse_pauling_data()
    if symbol in data:
        return data[symbol]
    _warn(f"Solid-state energy (Pauling regression) data for element {symbol} not found.")
    return None


@functools.cache
def _load_magpie_data() -> dict[str, dict]:
    data: dict[str, dict] = {}
    magpie_df = pd.read_csv(os.path.join(data_directory, "magpie.csv"))
    for _index, row in magpie_df.iterrows():
        key = row.iloc[0]
        data[key] = {
            "Number": int(row.iloc[1]),
            "MendeleevNumber": int(row.iloc[2]),
            "AtomicWeight": float(row.iloc[3]),
            "MeltingT": float(row.iloc[4]),
            "Column": int(row.iloc[5]),
            "Row": int(row.iloc[6]),
            "CovalentRadius": float(row.iloc[7]),
            "Electronegativity": float(row.iloc[8]),
            "NsValence": int(row.iloc[9]),
            "NpValence": int(row.iloc[10]),
            "NdValence": int(row.iloc[11]),
            "NfValence": int(row.iloc[12]),
            "NValence": int(row.iloc[13]),
            "NsUnfilled": int(row.iloc[14]),
            "NpUnfilled": int(row.iloc[15]),
            "NdUnfilled": int(row.iloc[16]),
            "NfUnfilled": int(row.iloc[17]),
            "NUnfilled": int(row.iloc[18]),
            "GSvolume_pa": float(row.iloc[19]),
            "GSbandgap": float(row.iloc[20]),
            "GSmagmom": float(row.iloc[21]),
            "SpaceGroupNumber": int(row.iloc[22]),
        }
    return data


def lookup_element_magpie_data(symbol: str, copy: bool = True):
    """
    Retrieve element data contained in the Magpie representation.

    Taken from Ward, L., Agrawal, A., Choudhary, A. et al.
    A general-purpose machine learning framework for
    predicting properties of inorganic materials.
    npj Comput Mater 2, 16028 (2016).
    https://doi.org/10.1038/npjcompumats.2016.28

    Args:
        symbol : the atomic symbol of the element to look up.
        copy: if True (default), return a copy of the data dictionary,
        rather than a reference to a cached object -- only use
        copy=False in performance-sensitive code and where you are
        certain the dictionary will not be modified!

    Returns:
        dict: Magpie feature dictionary for the element, or None if not found.

    """
    data = _load_magpie_data()
    if symbol in data:
        return data[symbol].copy() if copy else data[symbol]
    _warn(f"Magpie data for element {symbol} not found.")
    return None


@functools.cache
def _load_valence_data() -> dict[str, dict]:
    data: dict[str, dict] = {}
    valence_df = pd.read_csv(os.path.join(data_directory, "element_valence_modified.csv"))
    for _index, row in valence_df.iterrows():
        data[row.iloc[0]] = {"NValence": int(row.iloc[1])}
    return data


def lookup_element_valence_data(symbol: str, copy: bool = True):
    """
    Retrieve valence electron data.

    For d-block elements, the s and d electrons contribute to NValence.
    For p-block elements, the s and p electrons contribute to NValence.
    For s- and f-block elements, NValence is calculated from the Noble Gas electron configuration
        i.e.

    Args:
        symbol : the atomic symbol of the element to look up.
        copy: if True (default), return a copy of the data dictionary,
        rather than a reference to a cached object -- only use
        copy=False in performance-sensitive code and where you are
        certain the dictionary will not be modified!

    Returns:
        NValence (int): the number of valence electrons
        Returns None if the element was not found among the external
        data.
    """
    data = _load_valence_data()
    if symbol in data:
        return data[symbol].copy() if copy else data[symbol]
    _warn(f"Valence data for element {symbol} not found.")
    return None
