#!/usr/bin/env python

################################################################################
#  Copyright Keith T. Butler, Adam J. Jackson (2013)                           #
#                                                                              #
# This file is part of SMACT: lattice.py is free software: you can             #
# redistribute it and/or modify it under the terms of the GNU General Public   #
# License as published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.                          #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for    #
# more details.                                                                #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
################################################################################

from builtins import object
import numpy as np
import smact

class Lattice(object):
      """A unique set of Sites.

      Lattice objects define a general crystal structure, with a space group and
      a collection of Site objects. These Site objects have their own fractional
      coordinates and a list of possible oxidation states (see the Site class).

      Specific crystal structures with elements assigned to sites are
      "materials" and use the Atoms class from the Atomic Simulation
      Environment.

      Attributes:
          basis_sites: A list of Site objects [SiteA, SiteB, SiteC, ...]
          comprising the basis sites in Cartesian coordinates

          space_group: Integer space group number according to the
          International Tables for Crystallography.

          structurbericht:
          Structurbericht identity, if applicable (e.g. 'B1')

      Methods:
          lattice_vector_calc():

      """

      def __init__(self, sites, space_group=1, strukturbericht=False):
            self.sites = sites
            self.space_group = space_group
            self.strukturbericht = strukturbericht

class Site(object):
      """
      A single lattice site with a list of possible oxidation states.

      The Site object is primarily used within Lattice objects.

      Attributes:
          position: A list of fractional coordinates [x,y,z]
          oxidation_states: A list of possible oxidation states e.g. [-1,0,1]

      """

      def __init__(self, position, oxidation_states=[0]):
            self.position = position
            self.oxidation_states = oxidation_states

