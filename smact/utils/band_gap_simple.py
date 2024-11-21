"""Band gap simple."""

###############################################################################
# Copyright Daniel Davies, Adam J. Jackson (2013)                             #
#                                                                             #
# This file is part of SMACT: Band_gap_simple.py is free software: you can    #
# redistribute it and/or modify it under the terms of the GNU General Public  #
# License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version.                         #
# This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   #
# more details.                                                               #
# You should have received a copy of the GNU General Public License along with#
# this program.  If not, see <http://www.gnu.org/licenses/>.                  #
#                                                                             #
###############################################################################
from __future__ import annotations

from smact.properties import band_gap_Harrison

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compound band gap estimates from elemental data.")
    parser.add_argument("-a", "--anion", type=str, help="Element symbol for anion.")
    parser.add_argument("-c", "--cation", type=str, help="Element symbol for cation.")
    parser.add_argument("-d", "--distance", type=float, help="Internuclear separation.")
    parser.add_argument("-v", "--verbose", action="store_true", help="More Verbose output.")
    args = parser.parse_args()

    print(
        band_gap_Harrison(
            verbose=args.verbose,
            anion=args.anion,
            cation=args.cation,
            distance=args.distance,
        )
    )
