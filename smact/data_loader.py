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
import os

import pandas as pd

from smact import data_directory

# Module-level switch: print "verbose" warning messages
# about missing data.
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


def _get_data_rows(filename):
    """Generator for datafile entries by row for custom oxidation states lists. Skips rows with no oxidation states for performance."""
    with open(filename) as file:
        for line in file:
            line = line.strip()
            if line[0] != "#" and any(char.isdigit() for char in line):
                yield line.split()


def float_or_None(x):
    """Cast a string to a float or to a None."""
    try:
        return float(x)
    except ValueError:
        return None


_el_ox_states = None


def lookup_element_oxidation_states(symbol, copy=True):
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used is the SMACT default and
    most exhaustive list.

    Args:
    ----
        symbol (str) : the atomic symbol of the element to look up.
        copy (Optional(bool)): if True (default), return a copy of the
            oxidation-state list, rather than a reference to the cached
            data -- only use copy=False in performance-sensitive code
            and where the list will not be modified!

    Returns:
    -------
        list: List of known oxidation states for the element.

            Returns None if oxidation states for the Element were not
            found in the external data.

    """
    global _el_ox_states

    if _el_ox_states is None:
        _el_ox_states = {}
        for items in _get_data_rows(os.path.join(data_directory, "oxidation_states.txt")):
            _el_ox_states[items[0]] = [int(oxidationState) for oxidationState in items[1:]]

    if symbol in _el_ox_states:
        if copy:
            # _el_ox_states stores lists -> if copy is set, make an implicit
            # deep copy.  The elements of the lists are integers, which are
            # "value types" in Python.

            return list(_el_ox_states[symbol])
        else:
            return _el_ox_states[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Oxidation states for element {symbol} not found.")
        return None


_el_ox_states_icsd = None


def lookup_element_oxidation_states_icsd(symbol, copy=True):
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those found
    in the ICSD (and judged to be non-spurious).

    Args:
    ----
        symbol (str) : the atomic symbol of the element to look up.
        copy (Optional(bool)): if True (default), return a copy of the
            oxidation-state list, rather than a reference to the cached
            data -- only use copy=False in performance-sensitive code
            and where the list will not be modified!

    Returns:
    -------
        list: List of known oxidation states for the element.

            Return None if oxidation states for the Element were not
            found in the external data.

    """
    global _el_ox_states_icsd

    if _el_ox_states_icsd is None:
        _el_ox_states_icsd = {}

        for items in _get_data_rows(os.path.join(data_directory, "oxidation_states_icsd.txt")):
            _el_ox_states_icsd[items[0]] = [int(oxidationState) for oxidationState in items[1:]]
    if symbol in _el_ox_states_icsd:
        if copy:
            # _el_ox_states_icsd stores lists -> if copy is set, make an implicit
            # deep copy. The elements of the lists are integers, which are
            # "value types" in Python.
            return list(_el_ox_states_icsd[symbol])
        else:
            return _el_ox_states_icsd[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Oxidation states for element {symbol}not found.")
        return None


_el_ox_states_sp = None


def lookup_element_oxidation_states_sp(symbol, copy=True):
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those that
    are in the Pymatgen default lambda table for structure prediction.

    Args:
    ----
        symbol (str) : the atomic symbol of the element to look up.
        copy (Optional(bool)): if True (default), return a copy of the
            oxidation-state list, rather than a reference to the cached
            data -- only use copy=False in performance-sensitive code
            and where the list will not be modified!

    Returns:
    -------
        list: List of known oxidation states for the element.

            Return None if oxidation states for the Element were not
            found in the external data.

    """
    global _el_ox_states_sp

    if _el_ox_states_sp is None:
        _el_ox_states_sp = {}

        for items in _get_data_rows(os.path.join(data_directory, "oxidation_states_SP.txt")):
            _el_ox_states_sp[items[0]] = [int(oxidationState) for oxidationState in items[1:]]

    if symbol in _el_ox_states_sp:
        if copy:
            # _el_ox_states_sp stores lists -> if copy is set, make an implicit
            # deep copy.  The elements of the lists are integers, which are
            # "value types" in Python.

            return list(_el_ox_states_sp[symbol])
        else:
            return _el_ox_states_sp[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Oxidation states for element {symbol} not found.")
        return None


_el_ox_states_wiki = None


def lookup_element_oxidation_states_wiki(symbol, copy=True):
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those that
    are on Wikipedia (https://en.wikipedia.org/wiki/Template:List_of_oxidation_states_of_the_elements).

    Args:
    ----
        symbol (str) : the atomic symbol of the element to look up.
        copy (Optional(bool)): if True (default), return a copy of the
            oxidation-state list, rather than a reference to the cached
            data -- only use copy=False in performance-sensitive code
            and where the list will not be modified!

    Returns:
    -------
        list: List of known oxidation states for the element.

            Return None if oxidation states for the Element were not
            found in the external data.

    """
    global _el_ox_states_wiki

    if _el_ox_states_wiki is None:
        _el_ox_states_wiki = {}

        for items in _get_data_rows(os.path.join(data_directory, "oxidation_states_wiki.txt")):
            _el_ox_states_wiki[items[0]] = [int(oxidationState) for oxidationState in items[1:]]

    if symbol in _el_ox_states_wiki:
        if copy:
            # _el_ox_states_wiki stores lists -> if copy is set, make an implicit
            # deep copy.  The elements of the lists are integers, which are
            # "value types" in Python.

            return list(_el_ox_states_wiki[symbol])
        else:
            return _el_ox_states_wiki[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Oxidation states for element {symbol} not found.")
        return None


_el_ox_states_custom = None


def lookup_element_oxidation_states_custom(symbol, filepath, copy=True):
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list is specified by the user in a text file.

    Args:
    ----
        symbol (str) : the atomic symbol of the element to look up. "all" can be used to return the list (copy=True) or dict (copy=False) for all oxidation states.
        filepath (str) : the path to the text file containing the
            oxidation states data.
        copy (Optional(bool)): if True (default), return a copy of the
            oxidation-state list, rather than a reference to the cached
            data -- only use copy=False in performance-sensitive code
            and where the list will not be modified!

    Returns:
    -------
        list: List of known oxidation states for the element.

            Return None if oxidation states for the Element were not
            found in the external data.

    """
    global _el_ox_states_custom

    if _el_ox_states_custom is None:
        _el_ox_states_custom = {}
        for items in _get_data_rows(filepath):
            _el_ox_states_custom[items[0]] = [int(oxidationState) for oxidationState in items[1:]]
    if symbol in _el_ox_states_custom:
        if copy:
            # _el_ox_states_custom stores lists -> if copy is set, make an implicit
            # deep copy.  The elements of the lists are integers, which are
            # "value types" in Python.
            return list(_el_ox_states_custom[symbol])
        else:
            return _el_ox_states_custom[symbol]
    elif symbol == "all":
        if copy:
            return list(_el_ox_states_custom)
        else:
            return _el_ox_states_custom
    else:
        if _print_warnings:
            print(f"WARNING: Oxidation states for element {symbol} not found.")
        return None


_el_ox_states_icsd24 = None


def lookup_element_oxidation_states_icsd24(symbol, copy=True):
    """
    Retrieve a list of known oxidation states for an element.
    The oxidation states list used contains only those found
    in the 2024 version of the ICSD (and has >=5 reports).

    Args:
        symbol (str) : the atomic symbol of the element to look up.
        copy (Optional(bool)): if True (default), return a copy of the
            oxidation-state list, rather than a reference to the cached
            data -- only use copy=False in performance-sensitive code
            and where the list will not be modified!

    Returns:
        list: List of known oxidation states for the element.

            Returns None if oxidation states for the Element were not
            found in the external data.
    """
    global _el_ox_states_icsd24

    if _el_ox_states_icsd24 is None:
        _el_ox_states_icsd24 = {}

        for items in _get_data_rows(os.path.join(data_directory, "oxidation_states_icsd24_filtered.txt")):
            _el_ox_states_icsd24[items[0]] = [int(oxidationState) for oxidationState in items[1:]]

    if symbol in _el_ox_states_icsd24:
        if copy:
            # _el_ox_states_icsd24 stores lists -> if copy is set, make an implicit
            # deep copy.  The elements of the lists are integers, which are
            # "value types" in Python.

            return list(_el_ox_states_icsd24[symbol])
        else:
            return _el_ox_states_icsd24[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Oxidation states for element {symbol} not found.")
        return None


# Loader and cache for the element HHI scores.

_element_hhis = None


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
    global _element_hhis

    if _element_hhis is None:
        _element_hhis = {}

        with open(os.path.join(data_directory, "hhi.txt")) as file:
            for line in file:
                line = line.strip()

                if line[0] != "#":
                    items = line.split()

                    _element_hhis[items[0]] = (
                        float(items[1]),
                        float(items[2]),
                    )

    if symbol in _element_hhis:
        return _element_hhis[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: HHI data for element {symbol} not found.")

        return None


# Loader and cache for elemental data

_element_data = None


def lookup_element_data(symbol: str, copy: bool = True):
    """
    Retrieve tabulated data for an element.

    The table "data/element_data.txt" contains a collection of relevant
    atomic data. If a cache exists in the form of the module-level
    variable _element_data, this is returned. Otherwise, a dictionary is
    constructed from the data table and cached before returning it.

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
    global _element_data
    if _element_data is None:
        _element_data = {}
        keys = (
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
        for items in _get_data_rows(os.path.join(data_directory, "element_data.txt")):
            # First two columns are strings and should be left intact
            # Everything else is numerical and should be cast to a float
            # or, if not clearly a number, to None
            clean_items = items[0:2] + list(map(float_or_None, items[2:]))

            _element_data.update({items[0]: dict(list(zip(keys, clean_items, strict=False)))})

    if symbol in _element_data:
        if copy:
            # _element_open_babel_derived_data stores dictionaries
            # -> if copy is set, use the dict.copy() function to return
            # a copy. The values are all Python "value types", so
            # explicitly cloning the elements is not necessary to make
            # a deep copy.

            return _element_data[symbol].copy()
        else:
            return _element_data[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Elemental data for {symbol} not found.")
            print(_element_data)
        return None


# Loader and cache for the element Shannon radii datasets.

_element_shannon_radii_data = None


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
    global _element_shannon_radii_data

    if _element_shannon_radii_data is None:
        _element_shannon_radii_data = {}

        with open(os.path.join(data_directory, "shannon_radii.csv")) as file:
            reader = csv.reader(file)

            # Skip the first row (headers).

            next(reader)

            for row in reader:
                # For the shannon radii, there are multiple datasets for
                # different element/oxidation-state/coordination
                # combinations.

                key = row[0]

                dataset = {
                    "charge": int(row[1]),
                    "coordination": row[2],
                    "crystal_radius": float(row[3]),
                    "ionic_radius": float(row[4]),
                    "comment": row[5],
                }

                if key in _element_shannon_radii_data:
                    _element_shannon_radii_data[key].append(dataset)
                else:
                    _element_shannon_radii_data[key] = [dataset]

    if symbol in _element_shannon_radii_data:
        if copy:
            # _element_shannon_radii_data stores a list of dictionaries
            # -> if copy is set, copy the list and use the dict.copy()
            # function on each element.
            # The dictionary values are all Python "value types", so
            # nothing further is required to make a deep copy.
            return [item.copy() for item in _element_shannon_radii_data[symbol]]
        else:
            return _element_shannon_radii_data[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Shannon-radius data for element {symbol} not found.")

        return None


# Loader and cache for the machine-learned extended element Shannon radii datasets.

_element_shannon_radii_data_extendedML = None


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
        ionic_radius
            *float*
        comment
            *str*

    """
    global _element_shannon_radii_data_extendedML

    if _element_shannon_radii_data_extendedML is None:
        _element_shannon_radii_data_extendedML = {}

        with open(os.path.join(data_directory, "shannon_radii_ML_extended.csv")) as file:
            reader = csv.reader(file)

            # Skip the first row (headers).

            next(reader)

            for row in reader:
                # For the shannon radii, there are multiple datasets for
                # different element/oxidation-state/coordination
                # combinations.

                key = row[0]

                dataset = {
                    "charge": int(row[1]),
                    "coordination": row[2],
                    "crystal_radius": float(row[3]),
                    "ionic_radius": float(row[4]),
                    "comment": row[5],
                }

                if key in _element_shannon_radii_data_extendedML:
                    _element_shannon_radii_data_extendedML[key].append(dataset)
                else:
                    _element_shannon_radii_data_extendedML[key] = [dataset]

    if symbol in _element_shannon_radii_data_extendedML:
        if copy:
            # _element_shannon_radii_data_extendedML stores a list of dictionaries
            # -> if copy is set, copy the list and use the dict.copy()
            # function on each element.
            # The dictionary values are all Python "value types", so
            # nothing further is required to make a deep copy.
            return [item.copy() for item in _element_shannon_radii_data_extendedML[symbol]]
        else:
            return _element_shannon_radii_data_extendedML[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Extended Shannon-radius data for element {symbol} not found.")

        return None


# Loader and cache for the element solid-state energy (SSE) datasets.

_element_ssedata = None


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
    global _element_ssedata

    if _element_ssedata is None:
        _element_ssedata = {}

        with open(os.path.join(data_directory, "SSE.csv")) as file:
            reader = csv.reader(file)

            for row in reader:
                dataset = {
                    "AtomicNumber": int(row[1]),
                    "SolidStateEnergy": float(row[2]),
                    "IonisationPotential": float(row[3]),
                    "ElectronAffinity": float(row[4]),
                    "MullikenElectronegativity": float(row[5]),
                    "SolidStateRenormalisationEnergy": float(row[6]),
                }

                _element_ssedata[row[0]] = dataset

    if symbol in _element_ssedata:
        return _element_ssedata[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Solid-state energy data for element {symbol} not found.")

        return None


# Loader and cache for the revised (2015) element solid-state energy
# (SSE) datasets.

_element_sse2015_data = None


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
    global _element_sse2015_data

    if _element_sse2015_data is None:
        _element_sse2015_data = {}

        with open(os.path.join(data_directory, "SSE_2015.csv")) as file:
            reader = csv.reader(file)

            for row in reader:
                # Elements can have multiple SSE values depending on
                # their oxidation state

                key = row[0]

                dataset = {
                    "OxidationState": int(row[1]),
                    "SolidStateEnergy2015": float(row[2]),
                }

                if key in _element_sse2015_data:
                    _element_sse2015_data[key].append(dataset)
                else:
                    _element_sse2015_data[key] = [dataset]

    if symbol in _element_sse2015_data:
        if copy:
            return [item.copy() for item in _element_sse2015_data[symbol]]
        else:
            return _element_sse2015_data[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Solid-state energy (revised 2015) data for element {symbol} not found.")

        return None


# Loader and cache for the element solid-state energy (SSE) from Pauling
# electronegativity datasets.

_element_ssepauling_data = None


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
    global _element_ssepauling_data

    if _element_ssepauling_data is None:
        _element_ssepauling_data = {}

        with open(os.path.join(data_directory, "SSE_Pauling.csv")) as file:
            reader = csv.reader(file)

            for row in reader:
                dataset = {"SolidStateEnergyPauling": float(row[1])}

                _element_ssepauling_data[row[0]] = dataset

    if symbol in _element_ssepauling_data:
        return _element_ssepauling_data[symbol]
    else:
        if _print_warnings:
            print(
                "WARNING: Solid-state energy data from Pauling "
                " electronegativity regression fit for "
                f" element {symbol} not found."
            )

        return None


_element_magpie_data = None


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
        list:
            Magpie features.
        Returns None if the element was not found among the external
        data.

        Magpie features are dictionaries with the keys:




    """
    global _element_magpie_data

    if _element_magpie_data is None:
        _element_magpie_data = {}

        df = pd.read_csv(os.path.join(data_directory, "magpie.csv"))
        for _index, row in df.iterrows():
            key = row.iloc[0]

            dataset = {
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
            _element_magpie_data[key] = dataset

    if symbol in _element_magpie_data:
        return _element_magpie_data[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Magpie data for element {symbol} not found.")

        return None


_element_valence_data = None


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
    global _element_valence_data

    if _element_valence_data is None:
        _element_valence_data = {}

        df = pd.read_csv(os.path.join(data_directory, "element_valence_modified.csv"))
        for _index, row in df.iterrows():
            key = row.iloc[0]

            dataset = {"NValence": int(row.iloc[1])}
            _element_valence_data[key] = dataset

    if symbol in _element_valence_data:
        return _element_valence_data[symbol]
    else:
        if _print_warnings:
            print(f"WARNING: Valence data for element {symbol} not found.")

        return None
