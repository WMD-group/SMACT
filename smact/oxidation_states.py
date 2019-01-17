###############################################################################
# Copyright Daniel Davies (2019)                                              #
#                                                                             #
# This file is part of SMACT: oxidationstates.py is free software: you can    #
# redistribute it and/or modify it under the terms of the GNU General Public  #
# License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version.  This program is        #
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;   #
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A       #
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.   #
# You should have received a copy of the GNU General Public License along     #
# with this program.  If not, see <http://www.gnu.org/licenses/>.             #
#                                                                             #
###############################################################################

###############################################################################
# Collection of functions for the statistical analysis of oxidation states.   #
# It is possible to use the values obtained in the publication "Materials     #
# Discovery by Chemical Analogy: Role of Oxidation States in Structure        #
# Prediction" - DOI: 10.1039/C8FD00032H                                       #
# In future it will be possible to create a new probabilistic model using     #
# own database of materials.                                                  #
# and the functions in this module.                                           #
###############################################################################

from os import path
import json
from smact import data_directory, Element, Species

class Oxidationstate_probability_finder:
    '''


    '''

    def __init__(self, probability_table=None):
        '''
        Initialise the Oxidationstate_probability_finder class.

        Args:
            probability_table (dict): Lookup table to get probabilities for anion-cation pairs.
            Must be of the format {(anion,cation): probability, ...} e.g. {('F-1', 'Li1'): 1.0,...}.
            If none, the default table is loaded from the data directory.

        '''
        if probability_table == None:
            with open(path.join(data_directory,'oxidation_state_probability_table.json'), 'r') as f:
                probability_data = json.load(f)

            probability_table = {}
            for i in probability_data:
                probability_table[(i[0],i[1])] = i[2]

        self._probability_table = probability_table

        included_anions = set([i[0] for i in self._probability_table.keys()])
        included_cations = set([i[1] for i in self._probability_table.keys()])
        included_species = list(included_anions) + list(included_cations)

        self._included_species = included_species

    def get_pair_probability(self, species1, species2):
        '''
        Get the anion-cation oxidation state probability for a provided pair of smact Species.

        Args:
            species1 (smact.Species): Cation or anion species
            species2 (smact.Species): Cation or anion species

        '''

        # Check that there is one cation and one anion
        if (species1.oxidation > 0) and (species2.oxidation < 0):
            cation = species1
            anion = species2
        elif (species1.oxidation < 0) and (species2.oxidation > 0):
            anion = species1
            cation = species2
        else:
            raise ValueError("One cation and one anion required.")

        # Generate keys for lookup table
        cat_key = ''.join([cation.symbol,  str(cation.oxidation)])
        an_key = ''.join([anion.symbol, str(anion.oxidation)])

        # Check that both the species are included in the probability table
        if not all(elem in self._included_species for elem in [an_key, cat_key]):
            raise NameError("One or both of {0}, {1} are not in the probability table.".format(cat_key,an_key))


        prob = self._probability_table[(an_key,cat_key)]
        return prob

    def get_included_species(self):
        '''
        Returns a list of species for which there exists data in the probability table used.

        '''
        return self._included_species

