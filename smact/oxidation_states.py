"""
smact.oxidation_states: Module for predicting the likelihood of species
coexisting in a compound based on statistical analysis of oxidation states.
It is possible to use the values obtained in the publication Materials
Discovery by Chemical Analogy: Role of Oxidation States in Structure
Prediction - DOI: 10.1039/C8FD00032H.
"""

import json
from os import path
from typing import Dict, Optional, Tuple

from numpy import mean
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Specie as pmgSpecies

from smact import Element, Species, data_directory


class Oxidation_state_probability_finder:
    """
    Uses the model developed in the Faraday Discussions Paper (DOI:10.1039/C8FD00032H)
    to compute the likelihood of metal species existing in solids in the presence of certain anions.
    """

    def __init__(
        self, probability_table: Optional[Dict[Tuple[str, str], float]] = None
    ):
        """
        Args:
            probability_table (dict): Lookup table to get probabilities for anion-cation pairs.
                Must be of the format {(anion,cation): probability, ...} e.g. {('F-1', 'Li1'): 1.0,...}.
                If none, the default table is loaded from the data directory.
        """
        if probability_table == None:
            with open(
                path.join(
                    data_directory, "oxidation_state_probability_table.json"
                )
            ) as f:
                probability_data = json.load(f)
            # Put data into the required format
            probability_table = {}
            for i in probability_data:
                probability_table[(i[0], i[1])] = i[2]

        self._probability_table = probability_table
        # Define set of species for which we have data
        included_anions = {i[0] for i in self._probability_table.keys()}
        included_cations = {i[1] for i in self._probability_table.keys()}
        included_species = list(included_anions) + list(included_cations)

        self._included_species = included_species
        self._included_cations = included_cations
        self._included_anions = included_anions

    def _generate_lookup_key(self, species1: Species, species2: Species):
        """
        Internal function to generate keys to lookup table.

        Args:
            species1 (smact.Species): Species
            species2 (smact.Species): Species

        Returns:
            table_key (tuple): For looking up probability in the form (an_key, cat_key).

        """
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
        cat_key = "".join([cation.symbol, str(int(cation.oxidation))])
        an_key = "".join([anion.symbol, str(int(anion.oxidation))])

        # Check that both the species are included in the probability table
        if not all(
            elem in self._included_species for elem in [an_key, cat_key]
        ):
            raise NameError(
                f"One or both of [{cat_key}, {an_key}] are not in the probability table."
            )

        table_key = (an_key, cat_key)
        return table_key

    def pair_probability(self, species1: Species, species2: Species) -> float:
        """
        Get the anion-cation oxidation state probability for a provided pair of smact Species.
        i.e. :math:`P_{SA}=\\frac{N_{SX}}{N_{MX}}` in the original paper (DOI:10.1039/C8FD00032H).

        Args:
            species1 (smact.Species): Cation or anion species
            species2 (smact.Species): Cation or anion species

        Returns:
            prob (float): Species-anion probability

        """
        # Generate lookup table key and use it to look up probability
        probability_table_key = self._generate_lookup_key(species1, species2)
        prob = self._probability_table[probability_table_key]
        return prob

    def get_included_species(self):
        """
        Returns a list of species for which there exists data in the probability table used.
        """
        return self._included_species

    def compound_probability(
        self, structure: Structure, ignore_stoichiometry: bool = True
    ) -> float:
        """
        calculate overall probability for structure or composition.

        Args:
            structure (pymatgen.Structure): Compound for which the probability score will be generated.
                Can also be a list of pymatgen or SMACT Species.
            ignore_stoichiometry (bool): Whether to weight probabilities by stoichiometry.
                Defaults to false as decribed in the original paper.

        Returns:
            compound_prob (float): Compound probability
        """

        # Convert input to list of SMACT Species
        if type(structure) == list:
            if all(isinstance(i, Species) for i in structure):
                pass
            elif all(isinstance(i, pmgSpecies) for i in structure):
                structure = [Species(i.symbol, i.oxi_state) for i in structure]
            else:
                raise TypeError(
                    "Input requires a list of SMACT or Pymatgen species."
                )
        elif type(structure) == Structure:
            species = structure.species
            if not all(isinstance(i, pmgSpecies) for i in species):
                raise TypeError("Structure must have oxidation states.")
            else:
                structure = [Species(i.symbol, i.oxi_state) for i in structure]
        else:
            raise TypeError(
                "Input requires a list of SMACT or Pymatgen Species or a Structure."
            )

        # Put most electonegative element last in list by sorting by electroneg
        structure.sort(key=lambda x: x.pauling_eneg)

        # Define necessary species pairs
        anion = structure[-1]
        cations = [i for i in structure if i.oxidation > 0]
        species_pairs = [(anion, cation) for cation in cations]

        # Reduce down to unique pairs if ignoring stoichiometry
        if ignore_stoichiometry:
            species_pairs = list(set(species_pairs))

        # Do the maths
        pair_probs = [
            self.pair_probability(pair[0], pair[1]) for pair in species_pairs
        ]
        compound_prob = mean(pair_probs)
        return compound_prob
