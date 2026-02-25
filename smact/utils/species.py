"""Pure-string utilities for parsing and formatting SMACT species strings."""

from __future__ import annotations

import re

__all__ = ["get_sign", "parse_spec", "unparse_spec"]


def parse_spec(species: str) -> tuple[str, int]:
    """
    Parse a species string into its atomic symbol and oxidation state.

    :param species: the species string
    :return: a tuple of the atomic symbol and oxidation state

    """
    try:
        ele, oxi_state = re.match(r"([A-Za-z]+)([0-9]*[\+\-])", species).groups()
        # The regex guarantees oxi_state ends in '+' or '-', so this branch
        # is always taken; the else clause is unreachable by construction.
        charge = (int(oxi_state[:-1] or 1)) * (-1 if "-" in oxi_state else 1)
        return ele, charge
    except AttributeError:
        return _parse_spec_old(species)


def _parse_spec_old(species: str) -> tuple[str, int]:
    """
    Parse a species string into its atomic symbol and oxidation state.

    :param species: the species string
    :return: a tuple of the atomic symbol and oxidation state

    """
    match = re.match(r"[A-Za-z]+", species)
    if match is None:
        raise ValueError(f"Invalid species string (no element symbol found): {species!r}")
    ele = match.group(0)

    charge_match = re.search(r"\d+", species)
    ox_state = int(charge_match.group(0)) if charge_match else 0

    if "-" in species:
        ox_state *= -1

    # Handle cases of X+ or X- (instead of X1+ or X1-) as well as X0+ / X0-.
    # Short-circuit order matters: "0" wins over bare "+" or "-" so that "X0+"
    # correctly returns 0 rather than +1.
    if ox_state == 0 and "0" in species:
        ox_state = 0
    elif "+" in species and ox_state == 0:
        ox_state = 1
    elif ox_state == 0 and "-" in species:
        ox_state = -1

    return ele, ox_state


def unparse_spec(species: tuple[str, int], include_one: bool = True) -> str:
    """Unparse a species into a string representation.


    The analogue of :func:`parse_spec`.

    Args:
        species (tuple[str,int]): A tuple of (element, signed_charge).
        include_one (bool): If True, include charge of 1 in the output if charge is 1 or -1.


    Returns:
    -------
        String of {element}{absolute_charge}{sign}.

    Examples:
    --------
        >>> unparse_spec(("Fe", 2))
        'Fe2+'
        >>> unparse_spec(("O", -2))
        'O2-'

    """
    if include_one or abs(species[1]) != 1:
        return f"{species[0]}{abs(species[1])}{get_sign(species[1])}"
    else:
        return f"{species[0]}{get_sign(species[1])}"


def get_sign(charge: int) -> str:
    """
    Get string representation of a number's sign.

    Args:
    ----
        charge (int): The number whose sign to derive.

    Returns:
    -------
        sign (str): either '+', '-', or '' for neutral.

    """
    if charge > 0:
        return "+"
    elif charge < 0:
        return "-"
    else:
        return ""
