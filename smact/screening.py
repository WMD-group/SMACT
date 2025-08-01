"""A collection of tools for estimating physical properties
based on chemical composition.
"""

from __future__ import annotations

import itertools
import os
import warnings
from enum import Enum
from itertools import combinations
from typing import TYPE_CHECKING

from pymatgen.core import Composition

from smact import Element, element_dictionary, neutral_ratios
from smact.data_loader import (
    lookup_element_oxidation_states_custom as oxi_custom,
)
from smact.metallicity import metallicity_score
from smact.utils.composition import composition_dict_maker, formula_maker

try:
    from enum import StrEnum
except ImportError:

    class StrEnum(str, Enum):
        """Backport of Python 3.11's StrEnum for Python 3.10."""

        def __str__(self):
            return str(self.value)


if TYPE_CHECKING:
    import pymatgen


class SmactFilterOutputs(StrEnum):
    """Allowed outputs of the `smact_filter` function."""

    default = "default"
    formula = "formula"
    dict = "dict"


def pauling_test(
    oxidation_states: list[int],
    electronegativities: list[float],
    symbols: list[str] | None = None,
    repeat_anions: bool = True,
    repeat_cations: bool = True,
    threshold: float = 0.0,
):
    """
    Check if a combination of ions makes chemical sense,
        (i.e. positive ions should be of lower electronegativity).

    Args:
    ----
        oxidation_states (list):  oxidation states of elements in the compound
        electronegativities (list): the corresponding  Pauling electronegativities
            of the elements in the compound
        symbols (list) :  chemical symbols of each site
        threshold (float): a tolerance for the allowed deviation from
            the Pauling criterion
        repeat_anions : boolean, allow an anion to repeat in different
            oxidation states in the same compound
        repeat_cations : as above, but for cations

    Returns:
    -------
        bool:
            True if anions are more electronegative than
            cations, otherwise False

    """
    if symbols is None:
        symbols = []
    if repeat_anions and repeat_cations and threshold == 0.0:
        return eneg_states_test(oxidation_states, electronegativities)

    elif repeat_anions and repeat_cations:
        return eneg_states_test_threshold(oxidation_states, electronegativities, threshold=threshold)

    elif _no_repeats(
        oxidation_states,
        symbols,
        repeat_anions=repeat_anions,
        repeat_cations=repeat_cations,
    ):
        if threshold == 0.0:
            return eneg_states_test(oxidation_states, electronegativities)
        else:
            return eneg_states_test_threshold(oxidation_states, electronegativities, threshold=threshold)
    else:
        return False


def _no_repeats(
    oxidation_states: list[int],
    symbols: list[str],
    repeat_anions: bool = False,
    repeat_cations: bool = False,
):
    """
    Check if any anion or cation appears twice.

    Args:
    ----
        oxidation_states (list): oxidation states of species
        symbols (list): chemical symbols corresponding to oxidation
            states
        repeat_anions (bool): If True, anions may be repeated (e.g. O
            in -1 and -2 states)
        repeat_cations (bool): if True, cations may be repeated (e.g.
            Cu in +1 and +2 states)

    Returns:
    -------
        bool: True if no anion or cation is repeated, False otherwise
    """
    if repeat_anions is False and repeat_cations is False:
        return len(symbols) == len(set(symbols))
    else:
        anions, cations = [], []
        for state, symbol in zip(oxidation_states, symbols, strict=False):
            if state > 0:
                cations.append(symbol)
            else:
                anions.append(symbol)
        return not (
            (not repeat_anions and len(anions) != len(set(anions)))
            or (not repeat_cations and len(cations) != len(set(cations)))
        )


def pauling_test_old(
    ox: list[int],
    paul: list[float],
    symbols: list[str],
    repeat_anions: bool = True,
    repeat_cations: bool = True,
    threshold: float = 0.0,
):
    """
    Check if a combination of ions makes chemical sense,
    (i.e. positive ions should be of lower Pauling electronegativity).
    This function should give the same results as pauling_test but is
    not optimised for speed.

    Args:
    ----
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
    -------
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
    return not max(positive) - min(negative) > threshold


def eneg_states_test(ox_states: list[int], enegs: list[float]):
    """
    Internal function for checking electronegativity criterion.

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and eneg_states_test_alternate.

    Args:
    ----
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound

    Returns:
    -------
        bool : True if anions are more electronegative than
               cations, otherwise False

    """
    for (ox1, eneg1), (ox2, eneg2) in combinations(list(zip(ox_states, enegs, strict=False)), 2):
        if (
            eneg1 is None
            or eneg2 is None
            or ((ox1 > 0) and (ox2 < 0) and (eneg1 >= eneg2))
            or ((ox1 < 0) and (ox2 > 0) and (eneg1 <= eneg2))
        ):
            return False
    return True


def eneg_states_test_threshold(ox_states: list[int], enegs: list[float], threshold: float | None = 0):
    """
    Internal function for checking electronegativity criterion.

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and eneg_states_test_alternate.

    A 'threshold' option is added so that this constraint may be
    relaxed somewhat.

    Args:
    ----
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound
        threshold (Option(float)): a tolerance for the allowed deviation from
            the Pauling criterion

    Returns:
    -------
        bool : True if anions are more electronegative than
               cations, otherwise False

    """
    for (ox1, eneg1), (ox2, eneg2) in combinations(list(zip(ox_states, enegs, strict=False)), 2):
        if ((ox1 > 0) and (ox2 < 0) and ((eneg1 - eneg2) > threshold)) or (
            (ox1 < 0) and (ox2 > 0) and (eneg2 - eneg1) > threshold
        ):
            return False
    return True


def eneg_states_test_alternate(ox_states: list[int], enegs: list[float]):
    """
    Internal function for checking electronegativity criterion.

    This implementation appears to be slightly slower than
    eneg_states_test, but further testing is needed.

    Args:
    ----
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound

    Returns:
    -------
        bool : True if anions are more electronegative than
               cations, otherwise False

    """
    min_cation_eneg, max_anion_eneg = 10, 0
    for ox_state, eneg in zip(ox_states, enegs, strict=False):
        if ox_state < 1:
            min_cation_eneg = min(eneg, min_cation_eneg)
        else:
            max_anion_eneg = max(eneg, max_anion_eneg)
    return min_cation_eneg > max_anion_eneg


def ml_rep_generator(
    composition: list[Element] | list[str],
    stoichs: list[int] | None = None,
):
    """
    Function to take a composition of Elements and return a
    list of values between 0 and 1 that describes the composition,
    useful for machine learning.

    The list is of length 103 as there are 103 elements
    considered in total in SMACT.

    e.g. Li2O --> [0, 0, 2/3, 0, 0, 0, 0, 1/3, 0 ....]

    Inspired by the representation used by Legrain et al. DOI: 10.1021/acs.chemmater.7b00789

    Args:
    ----
        composition (list): Element objects in composition OR symbols of elements in composition
        stoichs (list): Corresponding stoichiometries in the composition

    Returns:
    -------
        norm (list): List of floats representing the composition that sum
            to one

    """
    if stoichs is None:
        stoichs = [1 for i, el in enumerate(composition)]

    ML_rep = [0 for i in range(1, 103)]
    if isinstance(composition[0], Element):
        for element, stoich in zip(composition, stoichs, strict=False):
            ML_rep[int(element.number) - 1] += stoich
    else:
        for element, stoich in zip(composition, stoichs, strict=False):
            ML_rep[int(Element(element).number) - 1] += stoich

    return [float(i) / sum(ML_rep) for i in ML_rep]


def smact_filter(
    els: tuple[Element] | list[Element],
    threshold: int | None = 8,
    stoichs: list[list[int]] | None = None,
    species_unique: bool = True,
    oxidation_states_set: str = "icsd24",
    return_output: SmactFilterOutputs = SmactFilterOutputs.default,
) -> list[tuple[str, int, int]] | list[tuple[str, int]] | list[str] | list[dict]:
    """Function that applies the charge neutrality and electronegativity
    tests in one go for simple application in external scripts that
    wish to apply the general 'smact test'.

     .. warning::
        For backwards compatibility in SMACT >=2.7, expllicitly set oxidation_states_set to 'smact14' if you wish to use the 2014 SMACT default oxidation states.
        In SMACT 3.0, the smact_filter function will be set to use a new default oxidation states set.

    Args:
    ----
        els (tuple/list): A list of smact.Element objects.
        threshold (int): Threshold for stoichiometry limit, default = 8.
        stoichs (list[int]): A selection of valid stoichiometric ratios for each site.
        species_unique (bool): Whether or not to consider elements in different oxidation states as unique in the results.
        oxidation_states_set (string): A string to choose which set of oxidation states should be chosen. Options are 'smact14', 'icsd16',"icsd24", 'pymatgen_sp' and 'wiki' for the  2014 SMACT default, 2016 ICSD, 2024 ICSD, pymatgen structure predictor and Wikipedia (https://en.wikipedia.org/wiki/Template:List_of_oxidation_states_of_the_elements) oxidation states respectively. A filepath to an oxidation states text file can also be supplied as well.
        return_output (SmactFilterOutputs): If set to 'default', the function will return a list of tuples containing the tuples of symbols, oxidation states and stoichiometry values. "Formula" returns a list of formulas and "dict" returns a list of dictionaries.

    Returns:
    -------
        allowed_comps (list): Allowed compositions for that chemical system
        in the form [(elements), (oxidation states), (ratios)] if species_unique=True and tuple=False
        or in the form [(elements), (ratios)] if species_unique=False and tuple=False.

    Example usage:
        >>> from smact.screening import smact_filter
        >>> from smact import Element
        >>> els = (Element("Cs"), Element("Pb"), Element("I"))
        >>> comps = smact_filter(els, threshold=5)
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
        >>> comps = smact_filter(els, stoichs=[[1], [1], [3]])
        >>> for comp in comps:
        >>>     print(comp)
        [('Cs', 'Pb', 'I'), (1, 2, -1), (1, 1, 3)]


    """
    # Get symbols and electronegativities
    symbols = tuple(e.symbol for e in els)
    electronegs = [e.pauling_eneg for e in els]

    # Select the specified oxidation states set:
    oxi_set = {
        "smact14": [e.oxidation_states_smact14 for e in els],
        "icsd16": [e.oxidation_states_icsd16 for e in els],
        "icsd24": [e.oxidation_states_icsd24 for e in els],
        "pymatgen_sp": [e.oxidation_states_sp for e in els],
        "wiki": [e.oxidation_states_wiki for e in els],
    }

    if oxidation_states_set in oxi_set:
        ox_combos = oxi_set[oxidation_states_set]
        if oxidation_states_set == "wiki":
            warnings.warn(
                "This set of oxidation states is sourced from Wikipedia. The results from using this set could be questionable and should not be used unless you know what you are doing and have inspected the oxidation states.",
                stacklevel=2,
            )
    elif os.path.exists(oxidation_states_set):
        ox_combos = (oxi_custom(e.symbol, oxidation_states_set) for e in els)
    else:
        raise (
            Exception(
                f'{oxidation_states_set} is not valid. Enter either "smact14", "icsd", "pymatgen","wiki" or a filepath to a textfile of oxidation states.'
            )
        )

    compositions = []
    for ox_states in itertools.product(*ox_combos):
        # Test for charge balance
        cn_e, cn_r = neutral_ratios(ox_states, stoichs=stoichs, threshold=threshold)
        # Electronegativity test
        if cn_e and pauling_test(ox_states, electronegs):
            for ratio in cn_r:
                compositions.append((symbols, ox_states, ratio))
    # Return list depending on whether we are interested in unique species combinations
    # or just unique element combinations.
    if species_unique:
        match return_output:
            case SmactFilterOutputs.default:
                return compositions
            case SmactFilterOutputs.formula:
                return [formula_maker(smact_filter_output=comp) for comp in compositions]
            case SmactFilterOutputs.dict:
                return [composition_dict_maker(smact_filter_output=comp) for comp in compositions]

    else:
        compositions = [(i[0], i[2]) for i in compositions]

        compositions = list(set(compositions))
        match return_output:
            case SmactFilterOutputs.default:
                return compositions
            case SmactFilterOutputs.formula:
                return [formula_maker(smact_filter_output=comp) for comp in compositions]
            case SmactFilterOutputs.dict:
                return [composition_dict_maker(smact_filter_output=comp) for comp in compositions]


# ---------------------------------------------------------------------
#                      Simplified SMACT Screening Logic
# ---------------------------------------------------------------------


def smact_validity(
    composition: pymatgen.core.Composition | str,
    use_pauling_test: bool = True,
    include_alloys: bool = True,
    check_metallicity: bool = False,
    metallicity_threshold: float = 0.7,
    oxidation_states_set: str | None = None,
    include_zero: bool = False,
    consensus: int = 3,
    commonality: str = "medium",
) -> bool:
    """
    Check if a composition is valid according to SMACT rules:
    1) Passes charge neutrality.
    2) Passes (optional) Pauling electronegativity test, or is considered an alloy or metal if so chosen.

    This function short-circuits, returning True as soon as a valid combination is found.

    Args:
        composition (Composition or str): Composition to check.
        use_pauling_test (bool): Whether to apply the Pauling EN test.
        include_alloys (bool): Consider pure metals valid automatically.
        check_metallicity (bool): If True, consider high metallicity valid.
        metallicity_threshold (float): Score threshold for metallicity validity.
        oxidation_states_set (str): Which set of oxidation states to use. If specified it overrides the making of the oxidation states set.
        include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.
        consensus (int): Minimum number of occurrences in literature for an ion to be considered valid. Default is 3.
        commonality (str): Excludes species below a certain proportion of appearances in literature with respect to the total number of reports of a given element (after the consensus threshold has been applied). "low" includes all species, "medium" excludes rare species below 10% occurrence, and "high" excludes non-majority species below 50% occurrence. "main" selects the species with the highest occurrence for a given element. Users may also specify their own threshold (float or int). Default is "medium".

    Returns:
        bool: True if the composition is valid, False otherwise.
    """
    from smact import _gcd_recursive, metals, neutral_ratios
    from smact.utils.oxidation import ICSD24OxStatesFilter

    if oxidation_states_set is not None and any([include_zero, consensus != 3, commonality != "medium"]):
        warnings.warn(
            "Parameters include_zero, consensus, and commonality are only used when oxidation_states_set is None",
            stacklevel=2,
        )

    if isinstance(composition, str):
        composition = Composition(composition)

    elem_symbols = tuple(composition.as_dict().keys())

    # Fast path for single elements
    if len(set(elem_symbols)) == 1:
        return True

    # Fast path for alloys
    if include_alloys and all(sym in metals for sym in elem_symbols):
        return True

    # Fast path for high metallicity compositions
    if check_metallicity:
        score = metallicity_score(composition)
        if score >= metallicity_threshold:
            return True

    # Convert composition counts -> stoichiometric ratios
    counts = [int(v) for v in composition.as_dict().values()]
    gcd_val = _gcd_recursive(*counts)
    stoichs = [(int(c // gcd_val),) for c in counts]
    threshold = max(int(c // gcd_val) for c in counts)

    # Build smact elements + electronegativities
    space = element_dictionary(elem_symbols)
    smact_elems = [e[1] for e in space.items()]
    electronegs = [e.pauling_eneg for e in smact_elems]

    # Get oxidation states data
    if oxidation_states_set is None:
        ox_filter = ICSD24OxStatesFilter()
        filtered_df = ox_filter.filter(consensus=consensus, include_zero=include_zero, commonality=commonality)
        oxidation_dict = {
            row["element"]: [int(x) for x in row["oxidation_state"].split()] for _, row in filtered_df.iterrows()
        }
        ox_combos = []
        for el in smact_elems:
            ox_el = oxidation_dict.get(el.symbol, None)
            if ox_el is not None:
                ox_combos.append(ox_el)
            else:
                return False
    elif oxidation_states_set == "smact14":
        ox_combos = [el.oxidation_states_smact14 for el in smact_elems]
    elif oxidation_states_set == "icsd16":
        ox_combos = [el.oxidation_states_icsd16 for el in smact_elems]
    elif oxidation_states_set == "pymatgen_sp":
        ox_combos = [el.oxidation_states_sp for el in smact_elems]
    elif oxidation_states_set == "icsd24" or oxidation_states_set is None:
        ox_combos = [el.oxidation_states_icsd24 for el in smact_elems]
    elif os.path.exists(oxidation_states_set):
        ox_combos = [oxi_custom(el.symbol, oxidation_states_set) for el in smact_elems]
    elif oxidation_states_set == "wiki":
        warnings.warn(
            "This set of oxidation states is from Wikipedia. Use with caution.",
            stacklevel=2,
        )
        ox_combos = [el.oxidation_states_wiki for el in smact_elems]
    else:
        raise ValueError(f"{oxidation_states_set} is not valid. Provide a known set or a valid file path.")

    # Check all possible oxidation state combinations
    for ox_states in itertools.product(*ox_combos):
        cn_e, cn_r = neutral_ratios(ox_states, stoichs=stoichs, threshold=threshold)

        if cn_e:
            if not use_pauling_test:
                return True

            try:
                en_ok = pauling_test(ox_states, electronegs)
            except TypeError:
                en_ok = True

            if en_ok:
                return True

    return False
