"""Miscellaneous tools for data parsing."""

import re
from typing import Tuple


def parse_spec(species: str) -> Tuple[str, int]:
    """Parse a species string into its element and charge."""
    ele = re.match(r"[A-Za-z]+", species).group(0)

    charge_match = re.search(r"\d+", species)
    charge = int(charge_match.group(0)) if charge_match else 0

    if '-' in species:
        charge *= -1

    return ele, charge


def unparse_spec(species: Tuple[str, int]) -> str:
    """Unparse a species into a string representation.

    The analogue of :func:`parse_spec`.

    """
    return f"{species[0]}{abs(species[1])}{get_sign(species[1])}"


def get_sign(charge: int) -> str:
    """Get string representation of a number's sign.

    Args:
        charge: The number whose sign to derive.

    Returns:
        Sign; either '+', '-' or '' for neutral.

    """
    if charge > 0:
        return '+'
    elif charge < 0:
        return '-'
    else:
        return ''
