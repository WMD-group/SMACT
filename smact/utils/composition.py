"""Utility functioms for handling elements, species, formulas and composition"""
import re
from collections import defaultdict


# Adapted from ElementEmbeddings and Pymatgen
def parse_formula(formula: str) -> dict[str, int]:
    """Parse a formula into a dict of el:amt

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
        expanded_sym = "".join(
            [f"{el}{amt}" for el, amt in unit_sym_dict.items()]
        )
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
