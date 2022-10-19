#!/usr/bin/env python
# Searches for site overlap between two lattices
# All input between lines 72 - 79
# Requires the surface_points.py file, which defines the different surfaces
################################################################################
# Copyright Keith T Butler    (2015)                                           #
#                                                                              #
# This file is part of SMACT: builder.py is free software: you can             #
# redistribute it and/or modify it under the terms of the GNU General Public   #
# License as published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.                          #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for    #
# more details.                                                                #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
#                                                                              #
################################################################################

from optparse import OptionParser

import numpy as np
import surface_points


def find_max_csl(surfs_1, surfs_2, multiplicity1, multiplicity2):
    """
    Given surface points and multiplicities of the surfaces this returns the maximal overlap fraction of the sites
    Attr:
    surfs : lists of the surface points on each side of the interface.
    multiplicity : lists of the multiplicity of the lattice vectors of u and v for each side of the interface.
    Returns:
    max_csl : float, the maximum fraction overlap found.
    """
    csl_values = []
    for surface_1 in surfs_1:
        if len(surface_1) > 0:
            surf_1_super = super_surface(
                np.asarray(surface_1), np.asarray(multiplicity1)
            )
            for surface_2 in surfs_2:
                if len(surface_2) > 0:
                    surf_2_super = super_surface(
                        np.asarray(surface_2), np.asarray(multiplicity2)
                    )
                    for i in np.arange(0, 1, 0.1):
                        for j in np.arange(0, 1, 0.1):
                            t_surf = translate(surf_2_super, [i, j])
                            csl_values.append(csl(surf_1_super, t_surf, multiplicity1))

    return max(csl_values)


def super_surface(surface, multiplicity):
    """Makes a super cell out of the surface coordinates"""

    surf_super = []
    for site in surface:
        for u in range(1, multiplicity[0] + 1):
            for v in range(1, multiplicity[1] + 1):
                surf_super.append(
                    [
                        (site[0] + (u - 1)) / multiplicity[0],
                        (site[1] + (v - 1)) / multiplicity[1],
                    ]
                )
    return np.asarray(surf_super)


def distance(a, b, mult):
    """Calculate separations, don't forget that we need to scale the separations by the multiplicity of the MAPI surface in each direction."""
    d1 = abs(a[0] - b[0])
    if d1 > 1:
        d1 = d1 - 1
    d2 = abs(a[1] - b[1])
    if d2 > 1:
        d2 = d2 - 1

    return np.sqrt((d1 * mult[0]) ** 2 + (d2 * mult[1]) ** 2)


def csl(surface1, surface2, mult_a, tol=0.15):
    """Takes two surfaces and calculates the number of co-inciding sites (within a tolerance)"""
    coincidence = 0.0
    for site_a in surface1:
        for site_b in surface2:
            if distance(site_a, site_b, mult_a) <= tol:
                coincidence = coincidence + 1.0
    return coincidence * 2 / (len(surface1) + len(surface2))


def wrapped(site):
    """Crude minimum image for this code"""
    if site[0] > 1:
        site[0] = site[0] - 1
    if site[1] > 1:
        site[1] = site[1] - 1
    if site[0] < 0:
        site[0] = site[0] + 1
    if site[1] < 0:
        site[1] = site[1] + 1
    return site


def translate(surface, T):
    """Translate the positions of the ions by a given vector"""
    for i, site in enumerate(surface):
        site = wrapped(site + T)
        surface[i] = site
    return surface


def get_comma_separated_args(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(","))


###### THESE ARE THE INPUT VARIABLES #######
parser = OptionParser()
parser.add_option(
    "-a",
    "--matera",
    action="store",
    type="string",
    dest="mater1",
    default="perovskite",
    help="The first material to consider",
)
parser.add_option(
    "-b",
    "--materb",
    action="store",
    type="string",
    dest="mater2",
    default="perovskite",
    help="The second material to consider",
)
parser.add_option(
    "-x",
    "--millera",
    action="store",
    type=str,
    dest="milla",
    default="100",
    help="The first materials miller index to consider, format : 100",
)
parser.add_option(
    "-y",
    "--millerb",
    action="store",
    type=str,
    dest="millb",
    default="100",
    help="The second materials miller index to consider, format : 100 ",
)
parser.add_option(
    "-u",
    "--multa",
    type="string",
    action="callback",
    dest="multa",
    callback=get_comma_separated_args,
    help="The first materials multiplicity, format : 2,2",
)
parser.add_option(
    "-v",
    "--multb",
    type="string",
    action="callback",
    dest="multb",
    callback=get_comma_separated_args,
    help="The second materials multiplicity, format : 3,3",
)


(options, args) = parser.parse_args()
material_1 = options.mater1
material_2 = options.mater2
surface_1 = str(options.milla)
surface_2 = str(options.millb)
multiplicity1 = list(map(int, options.multa))
multiplicity2 = list(map(int, options.multb))
######## INPUT OVER ###########
print(material_1, surface_1)
print(material_2, surface_2)
surfs_1 = getattr(surface_points, material_1)(surface_1)
surfs_2 = getattr(surface_points, material_2)(surface_2)

if len(surfs_1) > 0 and len(surfs_2) > 0:
    maxi = find_max_csl(surfs_1, surfs_2, multiplicity1, multiplicity2)
    print(maxi)
