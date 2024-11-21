"""Utility functioms for handling elements, species, formulas and composition."""

from __future__ import annotations

import re
from collections import defaultdict

from pymatgen.core import Composition

from smact.structure_prediction.utilities import unparse_spec


# Adapted from ElementEmbeddings and Pymatgen
def parse_formula(formula: str) -> dict[str, float]:
    """Parse a chemical formula into a dictionary of elements and their amounts.

    Args:
        formula (str): Chemical formula

    Returns:
        dict: Dictionary of element symbol: amount
    """
    regex = r"\(([^\(\)]+)\)\s*([\.e\d]*)"
    r = re.compile(regex)
    m = re.search(r, formula)
    if m:
        factor = 1.0
        if m.group(2) != "":
            factor = float(m.group(2))
        unit_sym_dict = _get_sym_dict(m.group(1), factor)
        expanded_sym = "".join([f"{el}{amt}" for el, amt in unit_sym_dict.items()])
        expanded_formula = formula.replace(m.group(), expanded_sym)
        return parse_formula(expanded_formula)
    return _get_sym_dict(formula, 1)


def _get_sym_dict(formula: str, factor: float) -> dict[str, float]:
    sym_dict: dict[str, float] = defaultdict(float)
    regex = r"([A-Z][a-z]*)\s*([-*\.e\d]*)"
    r = re.compile(regex)
    for m in re.finditer(r, formula):
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


def comp_maker(smact_filter_output: tuple[str, int, int] | tuple[str, int]) -> Composition:
    """Convert an item in the output of smact.screening.smact_filer into a Pymatgen Composition.

    Args:
        smact_filter_output (tuple[str, int, int]|tuple[str, int]): An item in the list returned from smact_filter

    Returns:
        composition (pymatgen.core.Composition): An instance of the Composition class
    """
    if len(smact_filter_output) == 2:
        form = []
        for el, ammt in zip(smact_filter_output[0], smact_filter_output[-1], strict=False):
            form.append(el)
            form.append(ammt)
        form = "".join(str(e) for e in form)
    else:
        form = {
            unparse_spec((el, ox)): ammt
            for el, ox, ammt in zip(
                smact_filter_output[0],
                smact_filter_output[1],
                smact_filter_output[2],
                strict=False,
            )
        }
    return Composition(form)


def formula_maker(smact_filter_output: tuple[str, int, int] | tuple[str, int]) -> str:
    """Convert an item in the output of smact.screening.smact_filter into a chemical formula.

    Args:
        smact_filter_output (tuple[str, int, int]|tuple[str, int]): An item in the list returned from smact_filter

    Returns:
            formula (str): A formula

    """
    return comp_maker(smact_filter_output).reduced_formula
