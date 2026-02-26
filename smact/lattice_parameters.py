"""
This module can be used to calculate roughly the lattice parameters of a
lattice type, based on the radii of the species on each site.
"""

from __future__ import annotations

import numpy as np

# Wurtzite geometry: tetrahedral bond angle ~ 109.47 deg = 2*arccos(-1/3)
_TETRAHEDRAL_HALF_ANGLE_SIN = 0.817  # sin(109.47 deg / 2)
_TETRAHEDRAL_COMPLEMENT_SIN = 0.335  # sin(109.47 deg - 90 deg)

# Body-centred tetragonal: empirical a/r ratio (β-Sn structure)
_BCT_A_FACTOR = 3.86

# Litharge (B10, PbO-type) empirical factors
_LITHARGE_SUM_FACTOR = 1.31
_LITHARGE_C_OVER_A = 1.26  # c/a from PbO data (http://www.mindat.org/min-2466.html)


def _check_positive(values: list[float] | float, name: str = "radius") -> None:
    """Raise ValueError if any radius value is not positive."""
    if isinstance(values, (int, float)):
        values = [values]
    if any(v <= 0 for v in values):
        raise ValueError(f"All {name} values must be positive, got {values}")


def cubic_perovskite(
    shannon_radius: list[float],
) -> tuple[float, float, float, float, float, float]:  # Cubic Pervoskite
    """
    The lattice parameters of the cubic perovskite structure.

    Args:
    ----
        shannon_radius (list) : The radii of the a,b,c ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(shannon_radius)
    limiting_factors = [2 * sum(shannon_radius[1:])]
    a = max(limiting_factors)
    b = a
    c = a
    #    space = a * np.sqrt(3) - 2 * shannon_radius[1]
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


def wurtzite(
    shannon_radius: list[float],
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of the wurtzite structure.

    Args:
    ----
        shannon_radius (list) : The radii of the a,b ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(shannon_radius)
    shannon_radius = sorted(shannon_radius, reverse=True)  # Geometry assumes atom A is larger
    # "Ideal" wurtzite structure
    # c/a = 1.633, u = 0.375
    #
    alpha = 90
    beta = 90
    gamma = 120
    #
    # Scenario A: A atoms are touching
    #   i.e. height is that of two tetrahegons with side length a
    #    = 2 * sqrt(2/3) * a
    if shannon_radius[0] > _TETRAHEDRAL_HALF_ANGLE_SIN * (shannon_radius[0] + shannon_radius[1]):
        a = 2 * shannon_radius[0]
        b = a
        c = 2 * np.sqrt(2.0 / 3.0) * a
    else:
        # Scenario B: regular wurtzite, similar sizes
        a = 2 * _TETRAHEDRAL_HALF_ANGLE_SIN * (shannon_radius[0] + shannon_radius[1])
        b = a
        c = (shannon_radius[0] + shannon_radius[1]) * (2 + 2 * _TETRAHEDRAL_COMPLEMENT_SIN)
    #    inner_space = a * (6**0.5) - (4*shannon_radius[0])
    return a, b, c, alpha, beta, gamma


# A1#
def fcc(
    covalent_radius: float,
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of the A1.

    Args:
    ----
        covalent_radius (list) : The radii of the a ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(covalent_radius)
    a = 2 * 2**0.5 * covalent_radius
    b = 2 * 2**0.5 * covalent_radius
    c = 2 * 2**0.5 * covalent_radius
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


# A2#
def bcc(
    covalent_radius: float,
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of the A2.

    Args:
    ----
        covalent_radius (list) : The radii of the a ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(covalent_radius)
    a = 4 * covalent_radius / np.sqrt(3)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


# A3#
def hcp(
    covalent_radius: float,
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of the hcp.

    Args:
    ----
        covalent_radius (list) : The radii of the a ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(covalent_radius)
    a = 2 * covalent_radius
    b = a
    c = (4.0 / 3.0) * 6**0.5 * covalent_radius
    alpha = 90
    beta = 90
    gamma = 120
    return a, b, c, alpha, beta, gamma


# A4#
def diamond(
    covalent_radius: float,
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of the diamond.

    Args:
    ----
        covalent_radius (list) : The radii of the a ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(covalent_radius)
    a = 8 * covalent_radius / np.sqrt(3)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


# A5#
def bct(
    covalent_radius: float,
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of the bct.

    Args:
    ----
        covalent_radius (list) : The radii of the a ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(covalent_radius)
    a = _BCT_A_FACTOR * covalent_radius
    b = a
    c = 2 * covalent_radius
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


# B1
def rocksalt(
    shannon_radius: list[float],
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of rocksalt.

    Args:
    ----
        shannon_radius (list) : The radii of the a,b ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(shannon_radius)
    limiting_factors = [
        2 * 2**0.2 * shannon_radius[0],
        2 * 2**0.2 * shannon_radius[1],
        2 * shannon_radius[0] + 2 * shannon_radius[1],
    ]
    a = max(limiting_factors)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


# B2
def b2(
    shannon_radius: list[float],
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of b2.

    Args:
    ----
        shannon_radius (list) : The radii of the a,b ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(shannon_radius)
    limiting_factors = [
        2 * (shannon_radius[0] + shannon_radius[0]) / np.sqrt(3),
        2 * shannon_radius[1],
        2 * shannon_radius[0],
    ]
    a = max(limiting_factors)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


# B3
def zincblende(
    shannon_radius: list[float],
) -> tuple[float, float, float, float, float, float]:
    """
    The lattice parameters of Zinc Blende.

    Args:
    ----
        shannon_radius (list) : The radii of the a,b ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(shannon_radius)
    limiting_factors = [
        2 * (max(shannon_radius) * np.sqrt(2)),
        4 * (shannon_radius[0] + shannon_radius[1]) ** (1.0 / 3.0),
    ]
    a = max(limiting_factors)
    b = a
    c = a
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


# Zn-S-Zn angle is ~109.5 degrees (from a tetrahedron). It is exactly 2*invCos(-1/3).
# The distance of that Zn-Zn (diagonally to half the face) is (using the cosine rule) is
# root[2(r1+r2)^2 - 2(r1+r2)^(2)cos(ZnSZn angle)].


# B10
def b10(
    shannon_radius: list[float],
) -> tuple[float, float, float, float, float, float]:  # Litharge
    """
    The lattice parameters of Litharge.

    Args:
    ----
        shannon_radius (list) : The radii of the a,b ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(shannon_radius)
    limiting_factors = [
        4 * (max(shannon_radius)) / np.sqrt(2),
        sum(shannon_radius) * _LITHARGE_SUM_FACTOR,
    ]
    a = max(limiting_factors)
    b = a
    c = a * _LITHARGE_C_OVER_A
    alpha = 90
    beta = 90
    gamma = 90
    return a, b, c, alpha, beta, gamma


def stuffed_wurtzite(
    shannon_radii: list[float],
) -> tuple[float, float, float, float, float, float]:
    """
    The stuffed wurtzite structure (e.g. LiGaGe) space group P63/mc.

    Args:
    ----
        shannon_radii (list) : The radii of the a,b,c ions

    Returns:
    -------
       (tuple):
           float values of lattics constants and
           angles (a, b, c, alpha, beta, gamma)

    """
    _check_positive(shannon_radii)
    rac = shannon_radii[2] + shannon_radii[1]
    x = rac * np.sin(np.radians(19.5))
    c = 2 * rac + x
    y = rac * np.sin(np.radians(70.5))
    a = y * np.sin(np.radians(120)) / np.sin(np.radians(30))
    b = a
    alpha = 90
    beta = 90
    gamma = 120
    return a, b, c, alpha, beta, gamma
