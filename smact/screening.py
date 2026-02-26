"""A collection of tools for estimating physical properties
based on chemical composition.
"""

from __future__ import annotations

import itertools
import os
import warnings
from itertools import combinations
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from collections.abc import Sequence

from pymatgen.core import Composition

from smact import Element, _gcd_recursive, element_dictionary, metals, neutral_ratios
from smact.data_loader import (
    lookup_element_oxidation_states_custom as oxi_custom,
)
from smact.metallicity import metallicity_score
from smact.utils.compat import StrEnum
from smact.utils.composition import composition_dict_maker, formula_maker
from smact.utils.oxidation import ICSD24OxStatesFilter

_NUM_ELEMENTS = 103

__all__ = [
    "SmactFilterOutputs",
    "eneg_states_test",
    "eneg_states_test_threshold",
    "ml_rep_generator",
    "pauling_test",
    "smact_filter",
    "smact_validity",
]

MIXED_VALENCE_ELEMENTS: frozenset[str] = frozenset(
    {
        # Transition metals
        "Fe",
        "Mn",
        "Co",
        "Cu",
        "Ni",
        "V",
        "Ti",
        "Cr",
        "Nb",
        "Mo",
        "W",
        "Re",
        "Ru",
        "Os",
        "Pd",
        "Ag",
        "Au",
        "Sn",
        "Sb",
        "Bi",
        # Lanthanides / actinides
        "Ce",
        "Eu",
        "Yb",
        "U",
    }
)


class SmactFilterOutputs(StrEnum):
    """Allowed outputs of the `smact_filter` function."""

    default = "default"
    formula = "formula"
    composition_dict = "composition_dict"


def _format_output(
    compositions: list,
    return_output: SmactFilterOutputs,
) -> list:
    """Format smact_filter compositions according to the requested output type.

    Args:
        compositions: List of composition tuples from smact_filter.
        return_output: The desired output format.

    Returns:
        Formatted list of compositions.
    """
    match return_output:
        case SmactFilterOutputs.default:
            return compositions
        case SmactFilterOutputs.formula:
            return [formula_maker(smact_filter_output=comp) for comp in compositions]
        case SmactFilterOutputs.composition_dict:
            return [composition_dict_maker(smact_filter_output=comp) for comp in compositions]
        case _:
            raise ValueError(f"Invalid return_output: {return_output}. Must be a SmactFilterOutputs value.")


_OXI_SET_ATTR_MAP: dict[str, str] = {
    "smact14": "oxidation_states_smact14",
    "icsd16": "oxidation_states_icsd16",
    "icsd24": "oxidation_states_icsd24",
    "pymatgen_sp": "oxidation_states_sp",
    "wiki": "oxidation_states_wiki",
}


def _get_oxidation_states(
    elements: Sequence[Element],
    oxidation_states_set: str,
) -> list[list[int] | None]:
    """Look up oxidation states for each element from the named set or a custom file.

    Args:
        elements: Sequence of smact.Element objects.
        oxidation_states_set: Name of a built-in set or a filepath to a custom file.

    Returns:
        List of oxidation state lists (may contain None for missing elements).

    Raises:
        ValueError: If the oxidation_states_set is not recognised and is not a valid filepath.
    """
    if oxidation_states_set in _OXI_SET_ATTR_MAP:
        attr = _OXI_SET_ATTR_MAP[oxidation_states_set]
        if oxidation_states_set == "wiki":
            warnings.warn(
                "This set of oxidation states is sourced from Wikipedia. The results from using this set could be "
                "questionable and should not be used unless you know what you are doing and have inspected the "
                "oxidation states.",
                stacklevel=3,
            )
        return [getattr(e, attr) for e in elements]
    if os.path.exists(oxidation_states_set):
        return cast("list[list[int] | None]", [oxi_custom(e.symbol, oxidation_states_set) for e in elements])
    raise ValueError(
        f'{oxidation_states_set} is not valid. Enter either "smact14", "icsd16", "icsd24", '
        '"pymatgen_sp", "wiki" or a filepath to a textfile of oxidation states.'
    )


def pauling_test(
    oxidation_states: Sequence[int],
    electronegativities: Sequence[float | None],
    symbols: Sequence[str] | None = None,
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
    oxidation_states: Sequence[int],
    symbols: Sequence[str],
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



def eneg_states_test(ox_states: Sequence[int], enegs: Sequence[float | None]):
    """
    Internal function for checking electronegativity criterion.

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and alternatives.

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


def eneg_states_test_threshold(ox_states: Sequence[int], enegs: Sequence[float | None], threshold: float = 0):
    """
    Internal function for checking electronegativity criterion.

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and alternatives.

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
        if eneg1 is None or eneg2 is None:
            return False
        if ((ox1 > 0) and (ox2 < 0) and ((eneg1 - eneg2) > threshold)) or (
            (ox1 < 0) and (ox2 > 0) and (eneg2 - eneg1) > threshold
        ):
            return False
    return True



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
        stoichs = [1] * len(composition)

    ML_rep = [0] * _NUM_ELEMENTS
    if isinstance(composition[0], Element):
        for element, stoich in zip(composition, stoichs, strict=False):
            ML_rep[int(element.number) - 1] += stoich  # type: ignore[union-attr]
    else:
        for element, stoich in zip(composition, stoichs, strict=False):
            ML_rep[int(Element(element).number) - 1] += stoich  # type: ignore[arg-type]

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
        return_output (SmactFilterOutputs): If set to 'default', the function will return a list of tuples containing the tuples of symbols, oxidation states and stoichiometry values. "formula" returns a list of formulas and "composition_dict" returns a list of dictionaries.

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
        ...     print(comp)
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
        ...     print(comp)
        [('Cs', 'Pb', 'I'), (1, 2, -1), (1, 1, 3)]


    """
    # Get symbols and electronegativities
    symbols = tuple(e.symbol for e in els)
    electronegs = [e.pauling_eneg for e in els]

    # Select the specified oxidation states set:
    ox_combos = _get_oxidation_states(els, oxidation_states_set)

    # Guard: raise early if any element has no oxidation states in the chosen set
    missing = [e.symbol for e, ox in zip(els, ox_combos, strict=False) if ox is None or len(ox) == 0]
    if missing:
        raise ValueError(
            f"No oxidation states found for {missing} in oxidation_states_set='{oxidation_states_set}'. "
            "Cannot enumerate charge-neutral compositions."
        )

    # After the guard, ox_combos contains only non-None, non-empty lists of ints
    ox_combos_typed = cast("list[list[int]]", ox_combos)

    compositions = []
    for ox_states in itertools.product(*ox_combos_typed):
        # Test for charge balance
        cn_r = neutral_ratios(list(ox_states), stoichs=stoichs, threshold=threshold)
        # Electronegativity test
        if cn_r and pauling_test(ox_states, electronegs):
            for ratio in cn_r:
                compositions.append((symbols, ox_states, ratio))
    # Return list depending on whether we are interested in unique species combinations
    # or just unique element combinations.
    if not species_unique:
        compositions = list({(i[0], i[2]) for i in compositions})
    return _format_output(compositions, return_output)


# ---------------------------------------------------------------------
#                      Simplified SMACT Screening Logic
# ---------------------------------------------------------------------


def smact_validity(
    composition: Composition | str,
    use_pauling_test: bool = True,
    include_alloys: bool = True,
    check_metallicity: bool = False,
    metallicity_threshold: float = 0.7,
    oxidation_states_set: str | None = None,
    include_zero: bool = False,
    consensus: int = 3,
    commonality: str = "medium",
    mixed_valence: bool = False,
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
        mixed_valence (bool): If True, allow mixed valence elements to be treated as separate species. Default is False.

    Returns:
        bool: True if the composition is valid, False otherwise.
    """
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
            row["element"]: [int(x) for x in str(row["oxidation_state"]).split()] for _, row in filtered_df.iterrows()
        }
        ox_combos = []
        for el in smact_elems:
            ox_el = oxidation_dict.get(el.symbol)
            if ox_el is not None:
                ox_combos.append(ox_el)
            else:
                return False
    else:
        ox_combos = _get_oxidation_states(smact_elems, oxidation_states_set)

    # Guard: if any element has no oxidation states in the chosen set (e.g. noble
    # gases in most databases), the composition cannot be charge-balanced.
    if any(ox is None or len(ox) == 0 for ox in ox_combos):
        return False

    # After the guard, ox_combos contains only non-None, non-empty lists of ints
    ox_combos_valid = cast("list[list[int]]", ox_combos)

    # Check all possible oxidation state combinations
    ox_valid = _is_valid_oxi_state(ox_combos_valid, stoichs, threshold, electronegs, use_pauling_test)

    if ox_valid:
        return True
    elif mixed_valence and any(el in MIXED_VALENCE_ELEMENTS for el in elem_symbols):
        # Guard against combinatorial blow-up before expanding
        projected = 1
        for el, ox, count in zip(elem_symbols, ox_combos_valid, stoichs, strict=False):
            projected *= len(ox) ** count[0] if el in MIXED_VALENCE_ELEMENTS else len(ox)
        if projected > 1_000_000:
            warnings.warn(
                "Mixed-valence expansion would generate too many combinations "
                f"({projected:,}); skipping to avoid excessive runtime.",
                stacklevel=2,
            )
            return False
        # Treat each atom of a mixed-valence element as an independent site.
        # threshold is computed from the original stoichs and remains valid after
        # expansion: expanded MV sites have stoich (1,) ≤ threshold, and
        # non-MV sites retain their original stoichs which are also ≤ threshold.
        ox_combos_valid, stoichs, electronegs = _expand_mixed_valence_comp(
            ox_combos_valid, stoichs, electronegs, elem_symbols
        )
        return _is_valid_oxi_state(ox_combos_valid, stoichs, threshold, electronegs, use_pauling_test)

    return False


def _expand_mixed_valence_comp(
    ox_combos: list[list[int]],
    stoichs: list[tuple[int, ...]],
    electronegs: list[float | None],
    elem_symbols: tuple[str, ...],
) -> tuple[list[list[int]], list[tuple[int, ...]], list[float | None]]:
    """Expand mixed-valence elements into individual single-stoichiometry sites."""
    new_ox_combos = []
    new_stoichs = []
    new_electronegs = []
    for el, ox, count, electroneg in zip(elem_symbols, ox_combos, stoichs, electronegs, strict=False):
        if el in MIXED_VALENCE_ELEMENTS:
            new_ox_combos.extend([ox] * count[0])
            new_electronegs.extend([electroneg] * count[0])
            new_stoichs.extend([(1,)] * count[0])
        else:
            new_ox_combos.append(ox)
            new_electronegs.append(electroneg)
            new_stoichs.append(count)
    return new_ox_combos, new_stoichs, new_electronegs


def _is_valid_oxi_state(
    ox_combos: list[list[int]],
    stoichs: Sequence[Sequence[int]],
    threshold: int,
    electronegs: list[float | None],
    use_pauling_test: bool = True,
) -> bool:
    """Return True if any oxidation-state combination satisfies charge neutrality and the Pauling criterion."""
    for ox_states in itertools.product(*ox_combos):
        cn_r = neutral_ratios(ox_states, stoichs=stoichs, threshold=threshold)

        if cn_r:
            if not use_pauling_test:
                return True

            try:
                en_ok = pauling_test(ox_states, electronegs)
            except TypeError:
                en_ok = True

            if en_ok:
                return True

    return False
