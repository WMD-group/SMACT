import itertools
import os
import warnings
from collections import namedtuple
from itertools import combinations
from typing import List, Optional, Tuple, Union

import numpy as np
import pymatgen
from pymatgen.core import Composition

import smact
from smact import Element, element_dictionary, neutral_ratios
from smact.data_loader import (
    lookup_element_oxidation_states_custom as oxi_custom,
)

# Use named tuple to improve readability of smact_filter outputs
_allowed_compositions = namedtuple(
    "Composition", ["element_symbols", "oxidation_states", "stoichiometries"]
)
_allowed_compositions_nonunique = namedtuple(
    "Composition", ["element_symbols", "stoichiometries"]
)


def pauling_test(
    oxidation_states: List[int],
    electronegativities: List[float],
    symbols: List[str] = [],
    repeat_anions: bool = True,
    repeat_cations: bool = True,
    threshold: float = 0.0,
):
    """Check if a combination of ions makes chemical sense,
        (i.e. positive ions should be of lower electronegativity).

    Args:
        ox (list):  oxidation states of elements in the compound
        paul (list): the corresponding  Pauling electronegativities
            of the elements in the compound
        symbols (list) :  chemical symbols of each site
        threshold (float): a tolerance for the allowed deviation from
            the Pauling criterion
        repeat_anions : boolean, allow an anion to repeat in different
            oxidation states in the same compound
        repeat_cations : as above, but for cations

    Returns:
        bool:
            True if anions are more electronegative than
            cations, otherwise False
    """

    if repeat_anions and repeat_cations and threshold == 0.0:
        return eneg_states_test(oxidation_states, electronegativities)

    elif repeat_anions and repeat_cations:
        return eneg_states_test_threshold(
            oxidation_states, electronegativities, threshold=threshold
        )

    else:
        if _no_repeats(
            oxidation_states,
            symbols,
            repeat_anions=repeat_anions,
            repeat_cations=repeat_cations,
        ):
            if threshold == 0.0:
                return eneg_states_test(oxidation_states, electronegativities)
            else:
                return eneg_states_test_threshold(
                    oxidation_states, electronegativities, threshold=threshold
                )
        else:
            return False


def _no_repeats(
    oxidation_states: List[int],
    symbols: List[str],
    repeat_anions: bool = False,
    repeat_cations: bool = False,
):
    """
    Check if any anion or cation appears twice.

    Args:
        oxidation_states (list): oxidation states of species
        symbols (list): chemical symbols corresponding to oxidation
            states
        repeat_anions (bool): If True, anions may be repeated (e.g. O
            in -1 and -2 states)
        repeat_cations (bool): if True, cations may be repeated (e.g.
            Cu in +1 and +2 states)

    Returns: bool
    """
    if repeat_anions is False and repeat_cations is False:
        return len(symbols) == len(set(symbols))
    else:
        anions, cations = [], []
        for state, symbol in zip(oxidation_states, symbols):
            if state > 0:
                cations.append(symbol)
            else:
                anions.append(symbol)
        if not repeat_anions and len(anions) != len(set(anions)):
            return False
        elif not repeat_cations and len(cations) != len(set(cations)):
            return False
        else:
            return True


def pauling_test_old(
    ox: List[int],
    paul: List[float],
    symbols: List[str],
    repeat_anions: bool = True,
    repeat_cations: bool = True,
    threshold: float = 0.0,
):
    """Check if a combination of ions makes chemical sense,
    (i.e. positive ions should be of lower Pauling electronegativity).
    This function should give the same results as pauling_test but is
    not optimised for speed.

    Args:
        ox (list):  oxidation states of the compound
        paul (list): the corresponding  Pauling electronegativities
            of the elements in the compound
        symbols (list) :  chemical symbols of each site.
        threshold (float): a tolerance for the allowed deviation from
            the Pauling criterion
        repeat_anions : boolean, allow an anion to repeat in different
            oxidation states in the same compound.
        repeat_cations : as above, but for cations.

    Returns:
        (bool):
            True if anions are more electronegative than
            cations, otherwise False

    """
    if None in paul:
        return False
    positive = []
    negative = []
    pos_ele = []
    neg_ele = []
    for i, state in enumerate(ox):
        if state > 0:
            if repeat_cations:
                positive.append(paul[i])
            elif symbols[i] in pos_ele:  # Reject materials where the same
                # cation occupies two sites.
                return False
            else:
                positive.append(paul[i])
                pos_ele.append(symbols[i])

        if state < 0:
            if repeat_anions:
                negative.append(paul[i])
            elif symbols[i] in neg_ele:  # Reject materials where the same
                # anion occupies two sites.
                return False
            else:
                negative.append(paul[i])
                neg_ele.append(symbols[i])

    if len(positive) == 0 or len(negative) == 0:
        return False
    if max(positive) == min(negative):
        return False
    if max(positive) - min(negative) > threshold:
        return False
    else:
        return True


def eneg_states_test(ox_states: List[int], enegs: List[float]):
    """
    Internal function for checking electronegativity criterion

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and eneg_states_test_alternate.

    Args:
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound

    Returns:
        bool : True if anions are more electronegative than
               cations, otherwise False

    """
    for (ox1, eneg1), (ox2, eneg2) in combinations(
        list(zip(ox_states, enegs)), 2
    ):
        if eneg1 is None or eneg2 is None:
            return False
        elif (ox1 > 0) and (ox2 < 0) and (eneg1 >= eneg2):
            return False
        elif (ox1 < 0) and (ox2 > 0) and (eneg1 <= eneg2):
            return False
    else:
        return True


def eneg_states_test_threshold(
    ox_states: List[int], enegs: List[float], threshold: Optional[float] = 0
):
    """Internal function for checking electronegativity criterion

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and eneg_states_test_alternate.

    A 'threshold' option is added so that this constraint may be
    relaxed somewhat.

    Args:
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound
        threshold (Option(float)): a tolerance for the allowed deviation from
            the Pauling criterion

    Returns:
        bool : True if anions are more electronegative than
               cations, otherwise False

    """
    for (ox1, eneg1), (ox2, eneg2) in combinations(
        list(zip(ox_states, enegs)), 2
    ):
        if (ox1 > 0) and (ox2 < 0) and ((eneg1 - eneg2) > threshold):
            return False
        elif (ox1 < 0) and (ox2 > 0) and (eneg2 - eneg1) > threshold:
            return False
    else:
        return True


def eneg_states_test_alternate(ox_states: List[int], enegs: List[float]):
    """Internal function for checking electronegativity criterion

    This implementation appears to be slightly slower than
    eneg_states_test, but further testing is needed.

    Args:
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound

    Returns:
        bool : True if anions are more electronegative than
               cations, otherwise False

    """
    min_cation_eneg, max_anion_eneg = 10, 0
    for ox_state, eneg in zip(ox_states, enegs):
        if ox_state < 1:
            min_cation_eneg = min(eneg, min_cation_eneg)
        else:
            max_anion_eneg = max(eneg, max_anion_eneg)
    return min_cation_eneg > max_anion_eneg


def ml_rep_generator(
    composition: Union[List[Element], List[str]],
    stoichs: Optional[List[int]] = None,
):
    """Function to take a composition of Elements and return a
    list of values between 0 and 1 that describes the composition,
    useful for machine learning.

    The list is of length 103 as there are 103 elements
    considered in total in SMACT.

    e.g. Li2O --> [0, 0, 2/3, 0, 0, 0, 0, 1/3, 0 ....]

    Inspired by the representation used by Legrain et al. DOI: 10.1021/acs.chemmater.7b00789

    Args:
        composition (list): Element objects in composition OR symbols of elements in composition
        stoichs (list): Corresponding stoichiometries in the composition

    Returns:
        norm (list): List of floats representing the composition that sum
            to one

    """
    if stoichs == None:
        stoichs = [1 for i, el in enumerate(composition)]

    ML_rep = [0 for i in range(1, 103)]
    if type(composition[0]) == Element:
        for element, stoich in zip(composition, stoichs):
            ML_rep[int(element.number) - 1] += stoich
    else:
        for element, stoich in zip(composition, stoichs):
            ML_rep[int(Element(element).number) - 1] += stoich

    norm = [float(i) / sum(ML_rep) for i in ML_rep]
    return norm


def smact_filter(
    els: Union[Tuple[Element], List[Element]],
    threshold: Optional[int] = 8,
    stoichs: Optional[List[List[int]]] = None,
    species_unique: bool = True,
    oxidation_states_set: str = "smact14",
    comp_tuple: bool = False,
) -> Union[List[Tuple[str, int, int]], List[Tuple[str, int]]]:
    """Function that applies the charge neutrality and electronegativity
    tests in one go for simple application in external scripts that
    wish to apply the general 'smact test'.

     .. warning::
        For backwards compatability in SMACT >=2.7, expllicitly set oxidation_states_set to 'smact14' if you wish to use the 2014 SMACT default oxidation states.
        In SMACT 3.0, the smact_filter function will be set to use a new default oxidation states set.

    Args:
        els (tuple/list): A list of smact.Element objects
        threshold (int): Threshold for stoichiometry limit, default = 8
        stoichs (list[int]): A selection of valid stoichiometric ratios for each site.
        species_unique (bool): Whether or not to consider elements in different oxidation states as unique in the results.
        oxidation_states_set (string): A string to choose which set of oxidation states should be chosen. Options are 'smact14', 'icsd', 'pymatgen' and 'wiki' for the  2014 SMACT default, 2016 ICSD, pymatgen structure predictor and Wikipedia (https://en.wikipedia.org/wiki/Template:List_of_oxidation_states_of_the_elements) oxidation states respectively. A filepath to an oxidation states text file can also be supplied as well.
        comp_tuple (bool): Whether or not to return the results as a named tuple of elements and stoichiometries (True) or as a normal tuple of elements and stoichiometries (False).
    Returns:
        allowed_comps (list): Allowed compositions for that chemical system
        in the form [(elements), (oxidation states), (ratios)] if species_unique=True
        or in the form [(elements), (ratios)] if species_unique=False.

    Example usage:
        >>> from smact.screening import smact_filter
        >>> from smact import Element
        >>> els = (Element('Cs'), Element('Pb'), Element('I'))
        >>> comps = smact_filter(els, threshold =5 )
        >>> for comp in comps:
        >>>     print(comp)
        [('Cs', 'Pb', 'I'), (1, -4, -1), (5, 1, 1)]
        [('Cs', 'Pb', 'I'), (1, 2, -1), (1, 1, 3)]
        [('Cs', 'Pb', 'I'), (1, 2, -1), (1, 2, 5)]
        [('Cs', 'Pb', 'I'), (1, 2, -1), (2, 1, 4)]
        [('Cs', 'Pb', 'I'), (1, 2, -1), (3, 1, 5)]
        [('Cs', 'Pb', 'I'), (1, 4, -1), (1, 1, 5)]

    Example (using stoichs):
        >>> from smact.screening import smact_filter
        >>> from smact import Element
        >>> comps = smact_filter(els, stoichs = [[1],[1],[3]], comp_tuple=True )
        >>> for comp in comps:
        >>>     print(comp)
        Composition(element_symbols=('Cs', 'Pb', 'I'), oxidation_states=(1, 2, -1), stoichiometries=(1, 1, 3))


    """

    compositions = []

    # Get symbols and electronegativities
    symbols = tuple(e.symbol for e in els)
    electronegs = [e.pauling_eneg for e in els]

    # Select the specified oxidation states set:
    oxi_set = {
        "smact14": [e.oxidation_states_smact14 for e in els],
        "icsd": [e.oxidation_states_icsd for e in els],
        "pymatgen": [e.oxidation_states_sp for e in els],
        "wiki": [e.oxidation_states_wiki for e in els],
    }
    if oxidation_states_set in oxi_set:
        ox_combos = oxi_set[oxidation_states_set]
    elif os.path.exists(oxidation_states_set):
        ox_combos = [oxi_custom(e.symbol, oxidation_states_set) for e in els]
    else:
        raise (
            Exception(
                f'{oxidation_states_set} is not valid. Enter either "smact14", "icsd", "pymatgen","wiki" or a filepath to a textfile of oxidation states.'
            )
        )
    if oxidation_states_set == "wiki":
        warnings.warn(
            "This set of oxidation states is sourced from Wikipedia. The results from using this set could be questionable and should not be used unless you know what you are doing and have inspected the oxidation states.",
            stacklevel=2,
        )

    for ox_states in itertools.product(*ox_combos):
        # Test for charge balance
        cn_e, cn_r = neutral_ratios(
            ox_states, stoichs=stoichs, threshold=threshold
        )
        # Electronegativity test
        if cn_e:
            electroneg_OK = pauling_test(ox_states, electronegs)
            if electroneg_OK:
                for ratio in cn_r:
                    compositions.append(
                        _allowed_compositions(symbols, ox_states, ratio)
                        if comp_tuple
                        else (symbols, ox_states, ratio)
                    )

    # Return list depending on whether we are interested in unique species combinations
    # or just unique element combinations.
    if species_unique:
        return compositions
    else:
        if comp_tuple:
            compositions = [
                _allowed_compositions_nonunique(i[0], i[2])
                for i in compositions
            ]
        else:
            compositions = [(i[0], i[2]) for i in compositions]
        compositions = list(set(compositions))
        return compositions


def smact_validity(
    composition: Union[pymatgen.core.Composition, str],
    use_pauling_test: bool = True,
    include_alloys: bool = True,
    oxidation_states_set: Union[str, bytes, os.PathLike] = "smact14",
) -> bool:
    """Check if a composition is valid according to the SMACT rules.

    Composition is considered valid if it passes the charge neutrality test and the Pauling electronegativity test.

     .. warning::
        For backwards compatability in SMACT >=2.7, expllicitly set oxidation_states_set to 'smact14' if you wish to use the 2014 SMACT default oxidation states.
        In SMACT 3.0, the smact_filter function will be set to use a new default oxidation states set.

    Args:
        composition (Union[pymatgen.core.Composition, str]): Composition/formula to check. This can be a pymatgen Composition object or a string.
        use_pauling_test (bool): Whether to use the Pauling electronegativity test
        include_alloys (bool): If True, compositions which only contain metal elements will be considered valid without further checks.
        oxidation_states_set (Union[str, bytes, os.PathLike]): A string to choose which set of
            oxidation states should be chosen for charge-balancing. Options are 'smact14', 'icsd',
            'pymatgen' and 'wiki' for the 2014 SMACT default, 2016 ICSD, pymatgen structure predictor and Wikipedia
            (https://en.wikipedia.org/wiki/Template:List_of_oxidation_states_of_the_elements) oxidation states respectively.
            A filepath to an oxidation states text file can also be supplied.

    Returns:
        bool: True if the composition is valid, False otherwise

    """
    if isinstance(composition, str):
        composition = Composition(composition)
    elem_symbols = tuple(composition.as_dict().keys())

    if len(set(elem_symbols)) == 1:
        return True
    if include_alloys:
        is_metal_list = [elem_s in smact.metals for elem_s in elem_symbols]
        if all(is_metal_list):
            return True

    count = tuple(composition.as_dict().values())
    count = [int(c) for c in count]
    # Reduce stoichiometry to gcd
    gcd = smact._gcd_recursive(*count)
    count = [int(c / gcd) for c in count]

    space = element_dictionary(elem_symbols)
    smact_elems = [e[1] for e in space.items()]
    electronegs = [e.pauling_eneg for e in smact_elems]

    if oxidation_states_set == "smact14" or oxidation_states_set is None:
        ox_combos = [e.oxidation_states_smact14 for e in smact_elems]
    elif oxidation_states_set == "icsd":
        ox_combos = [e.oxidation_states_icsd for e in smact_elems]
    elif oxidation_states_set == "pymatgen":
        ox_combos = [e.oxidation_states_sp for e in smact_elems]
    elif os.path.exists(oxidation_states_set):
        ox_combos = [
            oxi_custom(e.symbol, oxidation_states_set) for e in smact_elems
        ]
    elif oxidation_states_set == "wiki":
        warnings.warn(
            "This set of oxidation states is sourced from Wikipedia. The results from using this set could be questionable and should not be used unless you know what you are doing and have inspected the oxidation states.",
            stacklevel=2,
        )
        ox_combos = [e.oxidation_states_wiki for e in smact_elems]

    else:
        raise (
            Exception(
                f'{oxidation_states_set} is not valid. Enter either "smact14", "icsd", "pymatgen","wiki" or a filepath to a textfile of oxidation states.'
            )
        )

    threshold = np.max(count)
    compositions = []
    for ox_states in itertools.product(*ox_combos):
        stoichs = [(c,) for c in count]
        # Test for charge balance
        cn_e, cn_r = smact.neutral_ratios(
            ox_states, stoichs=stoichs, threshold=threshold
        )
        # Electronegativity test
        if cn_e:
            if use_pauling_test:
                try:
                    electroneg_OK = pauling_test(ox_states, electronegs)
                except TypeError:
                    # if no electronegativity data, assume it is okay
                    electroneg_OK = True
            else:
                electroneg_OK = True
            if electroneg_OK:
                for ratio in cn_r:
                    compositions.append(tuple([elem_symbols, ox_states, ratio]))
    compositions = [(i[0], i[2]) for i in compositions]
    compositions = list(set(compositions))
    if len(compositions) > 0:
        return True
    else:
        return False
