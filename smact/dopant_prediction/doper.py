from typing import Callable, List, Tuple, Type, Union

import numpy as np
from pymatgen.util import plotting

import smact
from smact import Element, element_dictionary
from smact.structure_prediction import mutation, utilities
from smact.structure_prediction.mutation import CationMutator


class Doper:
    """
    A class to search for n & p type dopants
    Methods: get_dopants, plot_dopants
    """

    def __init__(self, original_species: Tuple[str, ...], filepath: str = None):
        """
        Intialise the `Doper` class with a tuple of species

        Args:
            original_species: See :class:`~.Doper`.
            filepath (str): Path to a JSON file containing lambda table data.
        """
        self.original_species = original_species
        self.filepath = filepath
        self.results = None

    def _softmax(self, x):
        """
        Compute the softmax function for an input array.

        Args:
            x (numpy.ndarray): Input array.

        Returns:
            numpy.ndarray: Array of softmax values.
        """
        e_x = np.exp(x - np.max(x))
        return e_x / e_x.sum()

    def _apply_softmax(self, data: List[smact.Element]) -> List[smact.Element]:
        """
        Apply the softmax function to the probabilities in a list of Element objects.

        Args:
            data (List[smact.Element]): List of Element objects with probabilities.

        Returns:
            List[smact.Element]: List of Element objects with softmax-applied probabilities.
        """
        probabilities = [item[2] for item in data]
        softmax_probabilities = self._softmax(probabilities)

        for element, softmax_prob in zip(data, softmax_probabilities):
            element[2] = softmax_prob
        return data

    def _get_selectivity(
        self,
        data_list: List[smact.Element],
        cations: List[smact.Element],
        CM: Type[CationMutator],
        sub,
    ):
        data = data_list.copy()
        for dopants in data:
            if sub == "anion":
                dopants.append(1.0)
                continue
            selected_site, original_specie, sub_prob = dopants[:3]
            sum_prob = sub_prob
            for cation in cations:
                if cation != original_specie:
                    sum_prob += CM.sub_prob(cation, selected_site)

            selectivity = sub_prob / sum_prob
            selectivity = round(selectivity, 2)
            dopants.append(selectivity)

        return data

    def _get_dopants(
        self,
        element_objects: List[smact.Element],
        spicie_ions: List[str],
        ion_type: str,
    ):
        """
        Get possible dopants for a given list of elements and dopants.

        Args:
            element_objects (List[smact.Element]): List of Element objects.
            spicie_ions (List[str]): List of original species (anions or cations) as strings.
            ion_type (str): Identify which spicie to check.

        Returns:
            List[str]: List of possible dopants.
        """
        poss_n_type = set()
        poss_p_type = set()

        for element in element_objects:
            # i.e. element: "Zn", [-2, -1, 0, +1, +2]
            oxi_state = element.oxidation_states
            el_symbol = element.symbol
            for state in oxi_state:
                for ion in spicie_ions:
                    ele = utilities.unparse_spec((el_symbol, state))
                    _, charge = utilities.parse_spec(ion)

                    if ion_type == "anion":
                        if state > charge and state < 0:
                            poss_n_type.add(ele)
                        elif state < charge:
                            poss_p_type.add(ele)
                    elif ion_type == "cation":
                        if state > charge:
                            poss_n_type.add(ele)
                        elif state < charge and state > 0:
                            poss_p_type.add(ele)

        return list(poss_n_type), list(poss_p_type)

    def get_dopants(
        self, num_dopants: int = 5, apply_softmax=False, get_selectivity=True
    ) -> dict:
        """
        Args:
            num_dopants (int): The number of suggestions to return for n- and p-type dopants.
            apply_softmax (bool): Whether to apply softmax to probabilities. (default = True)
            get_selectivity (bool): Whether
        Returns:
            (dict): Dopant suggestions, given as a dictionary with keys
            "n_type_cation", "p_type_cation", "n_type_anion", "p_type_anion".

        Examples:
            >>> test = Doper(('Ti4+','O2-'))
            >>> print(test.get_dopants(num_dopants=2))
                {'n-type cation substitutions': [('Ta5+', 8.790371775858281e-05),
                ('Nb5+', 7.830035204694342e-05)],
                'p-type cation substitutions': [('Na1+', 0.00010060400812977031),
                ('Zn2+', 8.56373996146833e-05)],
                'n-type anion substitutions': [('F1-', 0.01508116810515677),
                ('Cl1-', 0.004737202729901607)],
                'p-type anion substitutions': [('N3-', 0.0014663800608945628),
                ('C4-', 9.31310255126729e-08)]}
        """

        cations, anions = [], []

        for ion in self.original_species:
            try:
                _, charge = utilities.parse_spec(ion)
                if charge > 0:
                    cations.append(ion)
                elif charge < 0:
                    anions.append(ion)
            except Exception as e:
                print(f"{e}: charge is not defined for {ion}!")

        CM = mutation.CationMutator.from_json(self.filepath)

        # call all elements
        element_objects = list(element_dictionary().values())

        poss_n_type_cat, poss_p_type_cat = self._get_dopants(
            element_objects, cations, "cation"
        )
        poss_n_type_an, poss_p_type_an = self._get_dopants(
            element_objects, anions, "anion"
        )

        n_type_cat, p_type_cat, n_type_an, p_type_an = [], [], [], []
        for cation in cations:
            cation_charge = utilities.parse_spec(cation)[1]

            for n_specie in poss_n_type_cat:
                n_specie_charge = utilities.parse_spec(n_specie)[1]
                if cation_charge >= n_specie_charge:
                    continue
                n_type_cat.append(
                    [n_specie, cation, CM.sub_prob(cation, n_specie)]
                )

            for p_specie in poss_p_type_cat:
                p_specie_charge = utilities.parse_spec(p_specie)[1]
                if cation_charge <= p_specie_charge:
                    continue
                p_type_cat.append(
                    [p_specie, cation, CM.sub_prob(cation, p_specie)]
                )

        for anion in anions:
            anion_charge = utilities.parse_spec(anion)[1]

            for n_specie in poss_n_type_an:
                n_specie_charge = utilities.parse_spec(n_specie)[1]
                if anion_charge >= n_specie_charge:
                    continue
                n_type_an.append(
                    [n_specie, anion, CM.sub_prob(anion, n_specie)]
                )

            for p_specie in poss_p_type_an:
                p_specie_charge = utilities.parse_spec(p_specie)[1]
                if anion_charge <= p_specie_charge:
                    continue
                p_type_an.append(
                    [p_specie, anion, CM.sub_prob(anion, p_specie)]
                )

        dopants_lists = [n_type_cat, p_type_cat, n_type_an, p_type_an]

        # sort by probability
        for dopants_list in dopants_lists:
            dopants_list.sort(key=lambda x: x[-1], reverse=True)

        # select top n elements
        dopants_lists = [
            dopants_list[:num_dopants] for dopants_list in dopants_lists
        ]

        if get_selectivity:
            for i in range(len(dopants_lists)):
                sub = "cation"
                if i > 1:
                    sub = "anion"
                dopants_lists[i] = self._get_selectivity(
                    dopants_lists[i], cations, CM, sub
                )

        if apply_softmax:
            # apply a softmax function
            for i in range(len(dopants_lists)):
                dopants_lists[i] = self._apply_softmax(dopants_lists[i])

        keys = [
            "n-type cation substitutions",
            "p-type cation substitutions",
            "n-type anion substitutions",
            "p-type anion substitutions",
        ]

        self.results = dict(zip(keys, dopants_lists))

        # return the top (num_dopants) results for each case
        return self.results

    def plot_dopants(self) -> None:
        """
        Plot the dopant suggestions using the periodic table heatmap.
        Args:
            None
        Returns:
            None
        """
        assert (
            self.results
        ), "Dopants are not calculated. Run get_dopants first."

        for dopant_type, dopants in self.results.items():
            dict_results = {
                utilities.parse_spec(x)[0]: y for x, _, y in dopants
            }
            plotting.periodic_table_heatmap(
                elemental_data=dict_results,
                cmap="rainbow",
                blank_color="gainsboro",
                edge_color="white",
            )
