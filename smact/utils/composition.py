"""Utility functions for handling elements, species, formulas and composition."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any

from pymatgen.core import Composition

from smact.structure_prediction.utilities import unparse_spec

_FORMULA_PAREN_RE = re.compile(r"\(([^\(\)]+)\)\s*([\.e\d]*)")
_SYMBOL_RE = re.compile(r"([A-Z][a-z]*)\s*([-*\.e\d]*)")

_MAX_FORMULA_DEPTH = 50
_STOICH_ONLY_TUPLE_LEN = 2


# Adapted from ElementEmbeddings and Pymatgen
def parse_formula(formula: str) -> dict[str, float]:
    """Parse a chemical formula into a dictionary of elements and their amounts.

    Uses iterative parenthesis expansion with a depth limit to prevent
    ``RecursionError`` on malformed input.

    Args:
        formula (str): Chemical formula

    Returns:
        dict: Dictionary of element symbol: amount

    Raises:
        ValueError: If the formula contains more than ``_MAX_FORMULA_DEPTH``
            levels of nested parentheses.
    """
    for _ in range(_MAX_FORMULA_DEPTH):
        m = _FORMULA_PAREN_RE.search(formula)
        if not m:
            break
        factor = float(m.group(2)) if m.group(2) != "" else 1.0
        unit_sym_dict = _get_sym_dict(m.group(1), factor)
        expanded_sym = "".join([f"{el}{amt}" for el, amt in unit_sym_dict.items()])
        formula = formula.replace(m.group(), expanded_sym)
    else:
        msg = f"Formula exceeds maximum nesting depth of {_MAX_FORMULA_DEPTH}: {formula!r}"
        raise ValueError(msg)
    return _get_sym_dict(formula, 1)


def _get_sym_dict(formula: str, factor: float) -> dict[str, float]:
    sym_dict: dict[str, float] = defaultdict(float)
    for m in _SYMBOL_RE.finditer(formula):
        el = m.group(1)
        amt = 1.0
        if m.group(2).strip() != "":
            amt = float(m.group(2))
        sym_dict[el] += amt * factor
        formula = formula.replace(m.group(), "", 1)
    if formula.strip():
        msg = f"{formula} is an invalid formula"
        raise ValueError(msg)

    return sym_dict


def comp_maker(smact_filter_output: tuple[Any, ...]) -> Composition:
    """Convert an item in the output of smact.screening.smact_filer into a Pymatgen Composition.

    Args:
        smact_filter_output: An item in the list returned from smact_filter.
            Either (elements, oxidation_states, stoichiometries) or (elements, stoichiometries).

    Returns:
        composition (pymatgen.core.Composition): An instance of the Composition class
    """
    if len(smact_filter_output) == _STOICH_ONLY_TUPLE_LEN:
        form_parts: list[str | int] = []
        for el, ammt in zip(smact_filter_output[0], smact_filter_output[-1], strict=True):
            form_parts.append(el)
            form_parts.append(ammt)
        form: str | dict[str, int] = "".join(str(e) for e in form_parts)
    else:
        form = {
            unparse_spec((el, ox)): ammt
            for el, ox, ammt in zip(
                smact_filter_output[0],
                smact_filter_output[1],
                smact_filter_output[2],
                strict=True,
            )
        }
    return Composition(form)


def formula_maker(smact_filter_output: tuple[Any, ...]) -> str:
    """Convert an item in the output of smact.screening.smact_filter into a chemical formula.

    Args:
        smact_filter_output: An item in the list returned from smact_filter.

    Returns:
            formula (str): A formula

    """
    return comp_maker(smact_filter_output).reduced_formula


def composition_dict_maker(smact_filter_output: tuple[Any, ...]) -> dict:
    """Convert an item in the output of smact.screening.smact_filter into a composition dictionary.

    Args:
        smact_filter_output: An item in the list returned from smact_filter.

    Returns:
        composition_dict (dict[str, float]): A composition dictionary
    """
    return comp_maker(smact_filter_output).as_dict()
