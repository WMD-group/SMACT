###############################################################################
# Copyright J. M. Skelton, D. W. Davies, A. J. Jackson (2016)                 #
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
Provide data from text files while transparently caching for efficiency

This module handles the loading of external data used to initialise the
core smact.Element and smact.Species classes.  It implements a
transparent data-caching system to avoid a large amount of I/O when
naively constructing several of these objects.  It also implements a
switchable system to print verbose warning messages about possible
missing data (mainly for debugging purposes).
"""

import csv
import os

from smact import data_directory

# Module-level switch: print "verbose" warning messages
# about missing data.
_print_warnings = False


def set_warnings(enable=True):
    """Set verbose warning messages on and off.

    In order to see any of the warnings, this function needs to be
    called _before_ the first call to the smact.Element()
    constructor.

    Args:
    enable (bool) : print verbose warning messages.
    """

    global _print_warnings
    _print_warnings = enable


def _get_data_rows(filename):
    """Generator for datafile entries by row"""
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line[0] != '#':
                yield line.split()

def float_or_None(x):
    """Cast a string to a float or to a None"""
    try:
        return float(x)
    except ValueError:
        return None

# Loader and cache for the element oxidation-state data.
_el_ox_states = None


def lookup_element_oxidation_states(symbol, copy=True):
    """
    Retrieve a list of known oxidation states for an element.

    Args:
        symbol (str) : the atomic symbol of the element to look up.
        copy (Optional(bool)): if True (default), return a copy of the
            oxidation-state list, rather than a reference to the cached
            data -- only use copy=False in performance-sensitive code
            and where the list will not be modified!

    Returns:
        list: List of known oxidation states for the element.

            Return None if oxidation states for the Element were not
            found in the external data.
    """

    global _el_ox_states

    if _el_ox_states is None:
        _el_ox_states = {}

        for items in _get_data_rows(os.path.join(data_directory,
                                                 "oxidation_states.txt")):
            _el_ox_states[items[0]] = [int(oxidationState)
                                       for oxidationState in items[1:]]

    if symbol in _el_ox_states:
        if copy:
            # _el_ox_states stores lists -> if copy is set, make an implicit
            # deep copy.  The elements of the lists are integers, which are
            # "value types" in Python.

            return [oxidationState for oxidationState in _el_ox_states[symbol]]
        else:
            return _el_ox_states[symbol]
    else:
        if _print_warnings:
            print("WARNING: Oxidation states for element {0} "
                  "not found.".format(symbol))
        return None

# Loader and cache for the element crustal abundances.

_element_crustal_abundances = None


def lookup_element_crustal_abundance(symbol):
    """
    Retrieve the crustal abundance for an element.

    Args:
        symbol (str) : the atomic symbol of the element to look up.

    Returns:
        The crustal abundance, or None if a value for the element was not found
        in the external data.
    """

    global _element_crustal_abundances

    if _element_crustal_abundances is None:
        _element_crustal_abundances = {}

        with open(os.path.join(data_directory, "crustal_abundance.txt"),
                  'r') as file:
            for line in file:
                line = line.strip()

                if line[0] != '#':
                    items = line.split()

                    _element_crustal_abundances[items[0]] = float(items[1])

    if symbol in _element_crustal_abundances:
        return _element_crustal_abundances[symbol]
    else:
        if _print_warnings:
            print("WARNING: Crustal-abundance data for element {0} "
                  "not found.".format(symbol))

        return None

# Loader and cache for the element HHI scores.

_element_hhis = None


def lookup_element_hhis(symbol):
    """
    Retrieve the HHI_R and HHI_p scores for an element.

    Args:
        symbol : the atomic symbol of the element to look up.

    Returns:
        A (HHI_p, HHI_R) tuple, or None if values for the elements were
        not found in the external data.
    """

    global _element_hhis

    if _element_hhis is None:
        _element_hhis = {}

        with open(os.path.join(data_directory, "HHIs.txt"),
                  'r') as file:
            for line in file:
                line = line.strip()

                if line[0] != '#':
                    items = line.split()

                    _element_hhis[items[0]] = (float(items[1]),
                                               float(items[2]))

    if symbol in _element_hhis:
        return _element_hhis[symbol]
    else:
        if _print_warnings:
            print("WARNING: HHI data for element "
                  "{0} not found.".format(symbol))

        return None

# Loader and cache for the Open Babel-derived element data.

_element_data = None

def lookup_element_data(symbol, copy=True):
    """
    Retrieve tabulated data for an element.

    The table "data/element_data.txt" contains a collection of relevant
    atomic data. If a cache exists in the form of the module-level
    variable _element_data, this is returned. Otherwise, a dictionary is
    constructed from the data table and cached before returning it.

    Args:
        symbol (str) : Atomic symbol for lookup

        copy (Optional(bool)) : if True (default), return a copy of the
            data dictionary, rather than a reference to the cached
            object -- only used copy=False in performance-sensitive code
            and where you are certain the dictionary will not be
            modified!

    Returns (dict): Dictionary of data for given element, keyed by
        column headings from data/element_data.txt
    """

    if _element_data is None:
        _element_data = {}
        keys = ('Symbol', 'Name', 'Z', 'Mass', 'r_cov', 'e_affinity',
                'p_eig', 's_eig', 'Abundance', 'el_neg', 'ion_pot')
        for items in _get_data_rows(os.path.join(data_directory,
                                                 "element.txt")):
            _element_data.update({items[1]: dict(zip(keys, items))})

_element_open_babel_derived_data = None

def lookup_element_open_babel_derived_data(symbol, copy=True):
    """
    Retrieve the Open Banel-derived data for an element.

    Args:
        symbol (str) : the atomic symbol of the element to look up.

    copy (Optional(bool)) : if True (default), return a copy of the data
        dictionary, rather than a reference to the cached object --
        only used copy=False in performance-sensitive code and where
        you are certain the dictionary will not be modified!

    Returns:
        A dictionary containing the data read from the Open Babel table,
        keyed with the column headings.
    """

    global _element_open_babel_derived_data

    if _element_open_babel_derived_data is None:
        _element_open_babel_derived_data = {}

        for items in _get_data_rows(os.path.join(data_directory,
                                                 "element.txt")):
            key = items[1]

            dataset = {}

            dataset['Number'] = int(items[0])

            are_neg = float(items[2])

            # these values are not (presently) used by SMACT
            # -> no need to print a warning.

            # if _print_warnings and are_neg == 0.0:
            #     print("WARNING: Aldred-Rochow electronegativity for "
            #           "element {0} may be set to the default value of"
            #           " zero in the data table.".format(symbol))

            dataset['ARENeg'] = are_neg

            r_cov = float(items[3])

            if _print_warnings and r_cov == 1.6:
                print("WARNING: Covalent radius for element {0} may "
                      "be set to the default value of 1.6 in the data"
                      " table.".format(key))

            dataset['RCov'] = float(items[3])

            dataset['RBO'] = float(items[4])

            r_vdw = float(items[5])

            # These values are not (presently) used by SMACT ->
            # no need to print a warning.

            # if _print_warnings and r_vdw == 2.0:
            #     print("WARNING: Van der Waals radius for element {0} "
            #           "may be set to the default value of 2.0 in the "
            #           " data table.".format(symbol))

            dataset['RVdW'] = float(r_vdw)

            max_bnd = int(items[6])

            # These values are not (presently) used by SMACT
            # -> no need to print a warning.

            # if _print_warnings and max_bnd == 6:
            #     print("WARNING: Maximum bond valence for element {0} "
            #           "may be set to the default value of 6 in the "
            #           "data table.".format(symbol))

            dataset['MaxBnd'] = max_bnd

            dataset['Mass'] = float(items[7])

            el_neg = float(items[8])

            if _print_warnings and el_neg == 0.0:
                print("WARNING: Pauling electronegativity for {0} may "
                      "be set to the default of zero in the "
                      "data table.".format(key))

            dataset['ElNeg.'] = el_neg

            ionization = float(items[9])

            if _print_warnings and ionization == 0.0:
                print("WARNING: Ionisation potential for {0} may be set"
                      " to the default of zero in the data "
                      "table.".format(key))

            dataset['Ionization'] = ionization

            el_affinity = float(items[10])

            if _print_warnings and el_affinity == 0.0:
                print("WARNING: Electron affinity for {0} may be set to"
                      " the default of zero in the "
                      "data table.".format(key))

            dataset['ElAffinity'] = el_affinity

            dataset['RGB'] = (float(items[11]), float(items[12]),
                              float(items[13]))
            dataset['Name'] = items[14]

            _element_open_babel_derived_data[items[1]] = dataset

    if symbol in _element_open_babel_derived_data:
        if copy:
            # _element_open_babel_derived_data stores dictionaries
            # -> if copy is set, use the dict.copy() function to return
            # a copy. The values are all Python "value types", so
            # explicitly cloning the elements is not necessary to make
            # a deep copy.

            return _element_open_babel_derived_data[symbol].copy()
        else:
            return _element_open_babel_derived_data[symbol]
    else:
        if _print_warnings:
            print("WARNING: Open Babel-derived element data for {0} not"
                  " found.".format(symbol))

        return None

# Loader and cache for the element eigenvalues.

_element_eigenvalues = None


def lookup_element_eigenvalue(symbol):
    """
    Retrieve the eigenvalue for an element.

    If the element is not found in the data table, or has a value of 0.0,
    a value of None is returned. The data table is "Eigenvalues.csv".

    Args:
        symbol (float): the atomic symbol of the element to look up.

    Returns:
        float: The eigenvalue, if a non-zero eigenvalue was found in the data.
            Otherwise, return None.
    """

    global _element_eigenvalues

    if _element_eigenvalues is None:
        _element_eigenvalues = {}

        with open(os.path.join(data_directory,
                               "Eigenvalues.csv"), 'r') as file:
            reader = csv.reader(file)

            # Skip the first row (headers).

            next(reader)

            for row in reader:
                _element_eigenvalues[row[0]] = float(row[1])

    if symbol in _element_eigenvalues:
        eig = _element_eigenvalues[symbol]
        if eig == 0.:
            return None
        else:
            return eig
    else:
        if _print_warnings:
            print("WARNING: Eigenvalue data for "
                  "element {0} not found.".format(symbol))

        return None

# Loader and cache for the element s eigenvalues.

_element_seigenvalues = None


def lookup_element_s_eigenvalue(symbol):
    """
    Retrieve the s eigenvalue for an element.

    If the element is not found in the data table, or has a value of 0.0,
    a value of None is returned. The data table is "Eigenvalues_s.csv".

    Args:
        symbol (float): the atomic symbol of the element to look up.

    Returns:
        float: s eigenvalue, if a non-zero eigenvalue was found in
            the external data. Otherwise, return None.
    """

    global _element_seigenvalues

    if _element_seigenvalues is None:
        _element_seigenvalues = {}

        with open(os.path.join(data_directory,
                               "Eigenvalues_s.csv"), 'rU') as file:
            reader = csv.reader(file)

            # Skip the first row (headers).

            next(reader)

            for row in reader:
                _element_seigenvalues[row[0]] = float(row[1])

    if symbol in _element_seigenvalues:
        eig = _element_seigenvalues[symbol]
        if eig == 0.:
            return None
        else:
            return eig

    else:
        if _print_warnings:
            print("WARNING: s-eigenvalue data for element {0}"
                  " not found.".format(symbol))

        return None

# Loader and cache for the element Shannon radii datasets.

_element_shannon_radii_data = None


def lookup_element_shannon_radius_data(symbol, copy=True):
    """
    Retrieve Shannon radii for known states of an element.

    Retrieve Shannon radii for known oxidation states and coordination
    environments of an element.

    Args:
        symbol (str) : the atomic symbol of the element to look up.

    copy (Optional(bool)): if True (default), return a copy of the data
        dictionary, rather than a reference to the cached object --
        only use copy=False in performance-sensitive code and where
        you are certain the dictionary will not be modified!

    Returns:
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

        with open(os.path.join(data_directory, "shannon_radii.csv"),
                  'rU') as file:
            reader = csv.reader(file)

            # Skip the first row (headers).

            next(reader)

            for row in reader:
                # For the shannon radii, there are multiple datasets for
                # different element/oxidation-state/coordination
                # combinations.

                key = row[0]

                dataset = {
                    'charge': int(row[1]),
                    'coordination': row[2],
                    'crystal_radius': float(row[3]),
                    'ionic_radius': float(row[4]),
                    'comment': row[5]
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
            return [item.copy() for item in
                    _element_shannon_radii_data[symbol]]
        else:
            return _element_shannon_radii_data[symbol]
    else:
        if _print_warnings:
            print("WARNING: Shannon-radius data for element {0} not "
                  "found.".format(symbol))

        return None

# Loader and cache for the element solid-state energy (SSE) datasets.

_element_ssedata = None


def lookup_element_sse_data(symbol):
    """
    Retrieve the solid-state energy (SSE) data for an element.

    Taken from J. Am. Chem. Soc., 2011, 133 (42), pp 16852-16960,
    DOI: 10.1021/ja204670s

    Args:
        symbol : the atomic symbol of the element to look up.

    Returns:
        A dictionary containing the SSE dataset for the element, or None
        if the element was not found among the external data.
    """

    global _element_ssedata

    if _element_ssedata is None:
        _element_ssedata = {}

        with open(os.path.join(data_directory,
                               "SSE.csv"), 'rU') as file:
            reader = csv.reader(file)

            for row in reader:
                dataset = {
                    'AtomicNumber': int(row[1]),
                    'SolidStateEnergy': float(row[2]),
                    'IonisationPotential': float(row[3]),
                    'ElectronAffinity': float(row[4]),
                    'MullikenElectronegativity': float(row[5]),
                    'SolidStateRenormalisationEnergy': float(row[6])
                    }

                _element_ssedata[row[0]] = dataset

    if symbol in _element_ssedata:
        return _element_ssedata[symbol]
    else:
        if _print_warnings:
            print("WARNING: Solid-state energy data for element {0} not"
                  " found.".format(symbol))

        return None

# Loader and cache for the revised (2015) element solid-state energy
# (SSE) datasets.

_element_sse2015_data = None


def lookup_element_sse2015_data(symbol, copy=True):
    """Retrieve SSE (2015) data for element in oxidation state.

    Retrieve the solid-state energy (SSE2015) data for an element in an
    oxidation state.  Taken from J. Solid State Chem., 2015, 231,
    pp138-144, DOI: 10.1016/j.jssc.2015.07.037

    Args:
    symbol : the atomic symbol of the element to look up.
    copy: if True (default), return a copy of the data dictionary,
        rather than a reference to a cached object -- only use
        copy=False in performance-sensitive code and where you are
        certain the dictionary will not be modified!

    Returns:
    A list of SSE datasets for the element, or None if the element was
    not found among the external data.
    """

    global _element_sse2015_data

    if _element_sse2015_data is None:
        _element_sse2015_data = {}

        with open(os.path.join(data_directory, "SSE_2015.csv"),
                  'rU') as file:
            reader = csv.reader(file)

            for row in reader:
                # Elements can have multiple SSE values depending on
                # their oxidation state

                key = row[0]

                dataset = {
                    'OxidationState': int(row[1]),
                    'SolidStateEnergy2015': float(row[2])}

                if key in _element_sse2015_data:
                    _element_sse2015_data[key].append(dataset)
                else:
                    _element_sse2015_data[key] = [dataset]

    if symbol in _element_sse2015_data:
        if copy:
            return [item.copy() for item in
                    _element_sse2015_data[symbol]]
        else:
            return _element_sse2015_data[symbol]
    else:
        if _print_warnings:
            print("WARNING: Solid-state energy (revised 2015) data for "
                  "element {0} not found.".format(symbol))

        return None

# Loader and cache for the element solid-state energy (SSE) from Pauling
# electronegativity datasets.

_element_ssepauling_data = None


def lookup_element_sse_pauling_data(symbol):
    """Retrieve Pauling SSE data

    Retrieve the solid-state energy (SSEPauling) data for an element
    from the regression fit when SSE2015 is plotted against Pauling
    electronegativity.  Taken from J. Solid State Chem., 2015, 231,
    pp138-144, DOI: 10.1016/j.jssc.2015.07.037

    Args:
    symbol (str) : the atomic symbol of the element to look up.

    Returns: A dictionary containing the SSE2015 dataset for the
        element, or None if the element was not found among the external
        data.
    """

    global _element_ssepauling_data

    if _element_ssepauling_data is None:
        _element_ssepauling_data = {}

        with open(os.path.join(data_directory, "SSE_Pauling.csv"),
                  'rU') as file:
            reader = csv.reader(file)

            for row in reader:
                dataset = {
                    'SolidStateEnergyPauling': float(row[1])}

                _element_ssepauling_data[row[0]] = dataset

    if symbol in _element_ssepauling_data:
        return _element_ssepauling_data[symbol]
    else:
        if _print_warnings:
            print("WARNING: Solid-state energy data from Pauling "
                  " electronegativity regression fit for "
                  " element {0} not found.".format(symbol))

        return None
