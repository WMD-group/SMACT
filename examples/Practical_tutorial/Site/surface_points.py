#!/usr/bin/env python
# Defines the positions (on a 2D grid) of the positions of under-coordinated
# atoms at the surface of a material
################################################################################
# Copyright Keith T Butler    (2015)		                               #
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

"""
The general form of the surface_points file.

Each new species or class of materials is defined as a function, the function contains a dictionary of the known surfaces of that material or class.
The surfaces contain the 2D projection of the under-coordinated surface atoms, there may be several possible terminations, then there are all listed.

The function has a list called 'exists', this contains all surfaces for which surface points are available. If you add a new surface be sure to include it in this list too.

Each dictionary item contains a list of the surface cuts and must end with an empty tuple, note the ,() at the end of each item.
"""


def anatase(miller):
    exists = ["001", "010", "110", "100", "101"]
    if miller in exists:
        surfaces = {}
        surfaces["110"] = (
            [0.0, 0.0],
            [0.2, 0],
            [0.3, 0.5],
            [0.5, 0.5],
            [0.7, 0.5],
            [0.8, 0.0],
        ), ()
        surfaces["001"] = ([0.5, 0.0], [0.5, 0.5]), ()
        surfaces["010"] = (
            [0.5, 0.0],
            [0.0, 0.0],
            [0.0, 0.2],
            [0.5, 0.5],
            [0.5, 0.75],
            [0.8, 0.0],
        ), ()
        surfaces["100"] = (
            [0.5, 0.0],
            [0.0, 0.0],
            [0.0, 0.2],
            [0.5, 0.5],
            [0.5, 0.75],
            [0.8, 0.0],
        ), ()
        surfaces["101"] = (
            (
                [0.0, 0.17],
                [0.5, 0.12],
                [0.0, 0.37],
                [0.0, 0.62],
                [0.5, 0.67],
                [0.5, 0.87],
            ),
            ([0.0, 0.14], [0.0, 0.34], [0.5, 0.68], [0.5, 0.84]),
            (
                [0.0, 0.15],
                [0.0, 0.36],
                [0.5, 0.4],
                [0.0, 0.57],
                [0.5, 0.65],
                [0.5, 0.86],
                [0.0, 0.9],
            ),
            (),
        )

        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def TiO2a(miller):
    exists = ["001", "010", "110", "100", "101"]
    if miller in exists:
        surfaces = {}
        surfaces["110"] = (
            [0.0, 0.0],
            [0.2, 0],
            [0.3, 0.5],
            [0.5, 0.5],
            [0.7, 0.5],
            [0.8, 0.0],
        ), ()
        surfaces["001"] = ([0.5, 0.0], [0.5, 0.5]), ()
        surfaces["010"] = (
            [0.5, 0.0],
            [0.0, 0.0],
            [0.0, 0.2],
            [0.5, 0.5],
            [0.5, 0.75],
            [0.8, 0.0],
        ), ()
        surfaces["100"] = (
            [0.5, 0.0],
            [0.0, 0.0],
            [0.0, 0.2],
            [0.5, 0.5],
            [0.5, 0.75],
            [0.8, 0.0],
        ), ()
        surfaces["101"] = (
            (
                [0.0, 0.17],
                [0.5, 0.12],
                [0.0, 0.37],
                [0.0, 0.62],
                [0.5, 0.67],
                [0.5, 0.87],
            ),
            ([0.0, 0.14], [0.0, 0.34], [0.5, 0.68], [0.5, 0.84]),
            (
                [0.0, 0.15],
                [0.0, 0.36],
                [0.5, 0.4],
                [0.0, 0.57],
                [0.5, 0.65],
                [0.5, 0.86],
                [0.0, 0.9],
            ),
            (),
        )
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def WO3(miller):
    exists = ["100", "110"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = (
            ([0.0, 0.25], [0, 0.75], [0.5, 0.25], [0.5, 0.75]),
            ([0.25, 0.2], [0.5, 0.25], [0.5, 0.75], [0.7, 0.7], [0.5, 1.0], [1.0, 1.0]),
            (),
        )
        surfaces["110"] = ([0.5, 0.3], [0.75, 0.3], [1.0, 0.8], [0.75, 0.80]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def perovskite(miller):
    exists = ["100", "110", "112"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0, 0], [0.5, 0.5]), ([0.5, 0], [0, 0.5], [0.5, 0.5]), ()
        surfaces["112"] = ([0.0, 0.0], [0.5, 0.0], [0.5, 0.5]), ()
        surfaces["110"] = ([0.0, 0.0], [0.0, 0.5], [0.75, 0.5]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def CH3NH3PbI3(miller):
    exists = ["100", "110", "112"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0, 0], [0.5, 0.5]), ([0.5, 0], [0, 0.5], [0.5, 0.5]), ()
        surfaces["112"] = ([0.0, 0.0], [0.5, 0.0], [0.5, 0.5]), ()
        surfaces["110"] = ([0.0, 0.0], [0.0, 0.5], [0.75, 0.5]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def SrTiO3(miller):
    exists = ["100", "110", "112"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0, 0], [0.5, 0.5]), ([0.5, 0], [0, 0.5], [0.5, 0.5]), ()
        surfaces["112"] = ([0.0, 0.0], [0.5, 0.0], [0.5, 0.5]), ()
        surfaces["110"] = ([0.0, 0.0], [0.0, 0.5], [0.75, 0.5]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def zincblende(miller):
    exists = ["100", "110"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0.75, 0.25], [0.0, 0.0]), ()
        surfaces["110"] = ([0.25, 0.9], [0.25, 0.4], [0.5, 0.7], [0.5, 0.2]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def CuIz(miller):
    exists = ["100", "110"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0.75, 0.25], [0.0, 0.0]), ()
        surfaces["001"] = ([0.75, 0.25], [0.0, 0.0]), ()
        surfaces["110"] = ([0.25, 0.9], [0.25, 0.4], [0.5, 0.7], [0.5, 0.2]), ()
        surfaces["011"] = ([0.25, 0.9], [0.25, 0.4], [0.5, 0.7], [0.5, 0.2]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def rocksalt(miller):
    exists = ["001", "100", "110", "011"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0.0, 0.0], [0.5, 0.5]), ()
        surfaces["001"] = ([0.0, 0.0], [0.5, 0.5]), ()
        surfaces["110"] = ([0.0, 0.0], [0.0, 0.5], [0.5, 0.0], [0.5, 0.5]), ()
        surfaces["011"] = ([0.0, 0.0], [0.0, 0.5], [0.5, 0.0], [0.5, 0.5]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def ZnTe(miller):
    exists = ["001", "100", "110", "011"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0.0, 0.0], [0.5, 0.5]), ()
        surfaces["001"] = ([0.0, 0.0], [0.5, 0.5]), ()
        surfaces["110"] = ([0.0, 0.0], [0.0, 0.5], [0.5, 0.0], [0.5, 0.5]), ()
        surfaces["011"] = ([0.0, 0.0], [0.0, 0.5], [0.5, 0.0], [0.5, 0.5]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def bixybite(miller):
    exists = ["100"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = (
            ([0.2, 0.9], [0.6, 0.9], [0.9, 0.6], [0.4, 0.4], [0.9, 0.4], [0.7, 0.1]),
            (
                [0.2, 0.2],
                [0.2, 0.7],
                [0.7, 0.2],
                [0.7, 0.7],
                [0.0, 0.3],
                [0.3, 0.5],
                [0.8, 0.5],
                [0.5, 0.8],
            ),
            (),
        )
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def rutile(miller):
    exists = ["001", "010", "110", "011"]
    if miller in exists:
        surfaces = {}
        surfaces["001"] = ([0.5, 0.5], [0.2, 0.8], [0.8, 0.2]), ()
        surfaces["010"] = ([0.5, 0.5], [0.7, 0.0]), ()
        surfaces["110"] = ([0.0, 0.9], [0.0, 0.45]), ()
        surfaces["011"] = ([0.0, 0.7], [0.3, 0.9], [0.2, 0.4], [0.5, 0.2]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def MoO3(miller):
    exists = ["001", "101"]
    if miller in exists:
        surfaces = {}
        surfaces["001"] = (
            ([0.9, 0.75], [0.7, 0.25], [0.6, 0.25], [0.4, 0.25]),
            (
                [0.9, 0.75],
                [0.7, 0.25],
                [0.6, 0.25],
                [0.6, 0.75],
                [0.4, 0.25],
                [0.4, 0.75],
                [0.3, 0.75],
                [0.0, 0.25],
            ),
            (),
        )
        surfaces["101"] = (
            [0.25, 1.0],
            [0.75, 0.7],
            [0.75, 0.66],
            [0.25, 0.5],
            [0.75, 0.3],
            [0.25, 0.1],
        ), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def wurtzite(miller):
    exists = ["100", "010", "110"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0, 0], [0, 0.37]), ()
        surfaces["010"] = ([0, 0], [0, 0.37]), ()
        surfaces["110"] = ([0, 0.8], [0.37, 0.8], [0.5, 0.17], [0.87, 0.17]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def GaN(miller):
    exists = ["100", "010", "110"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0, 0], [0, 0.37]), ()
        surfaces["010"] = ([0, 0], [0, 0.37]), ()
        surfaces["110"] = ([0, 0.8], [0.37, 0.8], [0.5, 0.17], [0.87, 0.17]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def SiC(miller):
    exists = ["100", "010", "110"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0, 0], [0, 0.37]), ()
        surfaces["010"] = ([0, 0], [0, 0.37]), ()
        surfaces["110"] = ([0, 0.8], [0.37, 0.8], [0.5, 0.17], [0.87, 0.17]), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def Cu2O(miller):
    exists = ["100", "110", "001", "011"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0.0, 0.0],), ()
        surfaces["110"] = ([0.0, 0.0],), ()
        surfaces["001"] = ([0.0, 0.0],), ()
        surfaces["011"] = ([0.0, 0.0],), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def In2S3(miller):
    exists = ["001", "010"]
    if miller in exists:
        surfaces = {}
        surfaces["001"] = (
            ([0.0, 0.25], [0.5, 0.0]),
            (
                [0.0, 0.0],
                [0.5, 0.0],
                [0.25, 0.75],
                [0.75, 0.75],
                [0.0, 0.5],
                [0.5, 0.0],
            ),
            ([0.0, 0.0], [0.75, 0.25], [0.0, 0.5], [0.25, 0.75], [0.25, 0.25]),
            (),
        )
        surfaces["010"] = (
            [0.5, 0.0],
            [0.75, 0.2],
            [0.25, 0.2],
            [0.75, 0.25],
            [0.75, 0.3],
            [0.25, 0.3],
            [0.0, 0.45],
            [0.75, 0.5],
            [0.25, 0.5],
            [0.75, 0.6, 0.25, 0.6],
            [0.75, 0.66],
            [0.25, 0.66],
            [0.5, 0.7],
            [0.0, 0.8],
            [0.25, 0.84],
            [0.75, 0.84],
            [0.75, 0.92],
            [0.75, 0.92],
        ), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def MnTiO3(miller):
    exists = ["010"]
    if miller in exists:
        surfaces = {}
        surfaces["010"] = (
            [0.6, 0.1],
            [0.9, 0.15],
            [0.2, 0.25],
            [0.9, 0.36],
            [0.52, 0.42],
            [0.25, 0.58],
            [0.9, 0.64],
            [0.6, 0.75],
            [0.9, 0.85],
            [0.23, 0.9],
        ), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def ZnTiO3(miller):
    exists = ["011"]
    if miller in exists:
        surfaces = {}
        surfaces["011"] = (
            [0.45, 0.09],
            [0.30, 0.45],
            [0.82, 0.26],
            [0.97, 0.60],
            [0.73, 0.82],
        ), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def SnS2(miller):
    exists = ["100"]
    if miller in exists:
        surfaces = {}
        surfaces["100"] = ([0.0, 0.0],), ([0.67, 0.33],), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def Ce2O3(miller):
    exists = ["101"]
    if miller in exists:
        surfaces = {}
        surfaces["101"] = ([0.69, 0.54],), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def LiNbO3(miller):
    exists = ["010"]
    if miller in exists:
        surfaces = {}
        surfaces["010"] = (
            [0.60, 0.90],
            [0.89, 0.83],
            [0.19, 0.73],
            [0.89, 0.61],
            [0.51, 0.57],
            [0.21, 0.40],
            [0.89, 0.33],
            [0.58, 0.23],
            [0.89, 0.11],
            [0.28, 0.07],
        ), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []


def Ce2S3(miller):
    exists = ["101", "110", "100", "010"]
    if miller in exists:
        surfaces = {}
        surfaces["101"] = (
            [0.25, 0.99],
            [0.25, 0.80],
            [0.25, 0.64],
            [0.25, 0.47],
            [0.75, 0.30],
            [0.25, 0.22],
            [0.75, 0.11],
        ), ()
        surfaces["110"] = (
            [0.93, 0.65],
            [0.89, 0.25],
            [0.78, 0.88],
            [0.70, 0.14],
            [0.72, 0.47],
            [0.54, 0.56],
            [0.43, 0.36],
            [0.39, 0.76],
            [0.29, 0.46],
            [0.20, 0.87],
            [0.04, 0.45],
        ), ()
        surfaces["100"] = (
            [0.75, 0.89],
            [0.75, 0.70],
            [0.25, 0.57],
            [0.75, 0.46],
            [0.75, 0.28],
            [0.25, 0.20],
            [0.75, 0.07],
        ), ()
        surfaces["010"] = (
            [0.73, 0.04],
            [0.13, 0.07],
            [0.85, 0.22],
            [0.35, 0.28],
            [0.23, 0.46],
            [0.63, 0.43],
            [0.99, 0.61],
            [0.36, 0.70],
            [0.86, 0.80],
            [0.49, 0.89],
        ), ()
        return surfaces[miller]
    else:
        print(
            "No non-polar surface",
            miller,
            "is currently in the database, maybe you want to add it.",
        )
    return []
