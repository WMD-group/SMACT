import os
from itertools import groupby
from typing import List, Optional, Tuple, Type

import numpy as np
from pymatgen.util import plotting
from tabulate import tabulate

import smact
from smact import data_directory
from smact.structure_prediction import mutation, utilities

SKIPSSPECIES_COSINE_SIM_PATH = os.path.join(
    data_directory,
    "species_rep/skipspecies_20221028_319ion_dim200_cosine_similarity.json",
)


class Doper:
    """
    A class to search for n & p type dopants
    Methods: get_dopants, plot_dopants
    """

    def __init__(
        self,
        original_species: Tuple[str, ...],
        filepath: Optional[str] = None,
        embedding: Optional[str] = None,
        use_probability: bool = True,
    ):
        """
        Initialise the `Doper` class with a tuple of species

        Args:
            original_species: See :class:`~.Doper`.
            filepath (str): Path to a JSON file containing lambda table data.
            embedding (str): Name of the species embedding to use. Currently only 'skipspecies' is supported.
            use_probability (bool): Whether to use the probability of substitution (calculated from `CationMutator`), or the raw similarity score/lambda value.

        """
        self.original_species = original_species
        self.filepath = filepath
        # filepath and embedding are mutually exclusive
        # check if both are provided
        if filepath and embedding:
            raise ValueError(
                "Only one of filepath or embedding should be provided"
            )
        if embedding and embedding != "skipspecies":
            raise ValueError(f"Embedding {embedding} is not supported")
        if embedding:
            self.cation_mutator = mutation.CationMutator.from_json(
                SKIPSSPECIES_COSINE_SIM_PATH
            )
        elif filepath:
            self.cation_mutator = mutation.CationMutator.from_json(filepath)
        else:
            # Default to Hautier data-mined lambda values
            self.cation_mutator = mutation.CationMutator.from_json(filepath)
        self.possible_species = list(self.cation_mutator.specs)
        self.lambda_threshold = self.cation_mutator.alpha("X", "Y")
        self.threshold = (
            1
            / self.cation_mutator.Z
            * np.exp(self.cation_mutator.alpha("X", "Y"))
        )
        self.use_probability = use_probability
        self.results = None

    def _get_selectivity(
        self,
        data_list: List[smact.Element],
        cations: List[smact.Element],
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
                    sum_prob += self._calculate_species_sim_prob(
                        cation, selected_site
                    )

            selectivity = sub_prob / sum_prob
            selectivity = round(selectivity, 2)
            dopants.append(selectivity)
            assert len(dopants) == 4
        return data

    def _merge_dicts(self, keys, dopants_list, groupby_list):
        merged_dict = dict()
        for k, dopants, group in zip(keys, dopants_list, groupby_list):
            merged_values = dict()
            merged_values["sorted"] = dopants
            for key, value in group.items():
                merged_values[key] = sorted(
                    value, key=lambda x: x[2], reverse=True
                )
            merged_dict[k] = merged_values
        return merged_dict

    def _get_dopants(
        self,
        specie_ions: List[str],
        ion_type: str,
    ):
        """
        Get possible dopants for a given list of elements and dopants.

        Args:
            specie_ions (List[str]): List of original species (anions or cations) as strings.
            ion_type (str): Identify which species to check.

        Returns:
            List[str]: List of possible dopants.
        """
        poss_n_type = set()
        poss_p_type = set()
        for spec in self.possible_species:
            _, state = utilities.parse_spec(spec)
            for ion in specie_ions:
                _, charge = utilities.parse_spec(ion)

                if ion_type == "anion":
                    if state > charge and state < 0:
                        poss_n_type.add(spec)
                    elif state < charge:
                        poss_p_type.add(spec)
                elif ion_type == "cation":
                    if state > charge:
                        poss_n_type.add(spec)
                    elif state < charge and state > 0:
                        poss_p_type.add(spec)

        return list(poss_n_type), list(poss_p_type)

    def _calculate_species_sim_prob(self, species1, species2):
        """
        Calculate the similarity/probability between two species.

        Args:
            species1 (str): The first species.
            species2 (str): The second species.

        Returns:
            float: The similarity between the two species.
        """
        if self.use_probability:
            return self.cation_mutator.sub_prob(species1, species2)
        else:
            return self.cation_mutator.get_lambda(species1, species2)

    def get_dopants(
        self, num_dopants: int = 5, get_selectivity=True, group_by_charge=True
    ) -> dict:
        """
        Args:
            num_dopants (int): The number of suggestions to return for n- and p-type dopants.
            get_selectivity (bool): Whether to calculate the selectivity of the dopants.
            group_by_charge (bool): Whether to group the dopants by charge.
        Returns:
            (dict): Dopant suggestions, given as a dictionary with keys
            "n-type cation substitutions", "p-type cation substitutions", "n-type anion substitutions", "p-type anion substitutions".

        Examples:
            >>> test = Doper(('Ti4+','O2-'))
            >>> print(test.get_dopants(num_dopants=2))
                {'n-type anion substitutions': {'-1': [['F1-', 'O2-', 0.01508116810515677, 1.0],
                                       ['Cl1-','O2-', 0.004737202729901607, 1.0]],
                                'sorted': [['F1-',
                                            'O2-',
                                            0.01508116810515677,
                                            1.0],
                                           ['Cl1-',
                                            'O2-',
                                            0.004737202729901607,
                                            1.0]]},
                'n-type cation substitutions': {'5': [['Ta5+',
                                        'Ti4+',
                                        8.790371775858281e-05,
                                        1.0],
                                       ['Nb5+',
                                        'Ti4+',
                                        7.830035204694342e-05,
                                        1.0]],
                                 '6': [['W6+',
                                        'Ti4+',
                                        3.4638026110457894e-05,
                                        1.0],
                                       ['Mo6+',
                                        'Ti4+',
                                        1.6924395455176864e-05,
                                        1.0]],
                                 'sorted': [['Ta5+',
                                             'Ti4+',
                                             8.790371775858281e-05,
                                             1.0],
                                            ['Nb5+',
                                             'Ti4+',
                                             7.830035204694342e-05,
                                             1.0]]},
                'p-type anion substitutions': {'-3': [['N3-',
                                        'O2-',
                                        0.0014663800608945628,
                                        1.0]],
                                'sorted': [['N3-',
                                            'O2-',
                                            0.0014663800608945628,
                                            1.0]]},
                'p-type cation substitutions': {'1': [['Na1+',
                                        'Ti4+',
                                        0.00010060400812977031,
                                        1.0],
                                       ['Li1+',
                                        'Ti4+',
                                        4.90559802023167e-05,
                                        1.0]],
                                 '2': [['Zn2+',
                                        'Ti4+',
                                        8.56373996146833e-05,
                                        1.0],
                                       ['Mn2+',
                                        'Ti4+',
                                        8.563568688381837e-05,
                                        1.0]],
                                 '3': [['Fe3+',
                                        'Ti4+',
                                        6.259479321178562e-05,
                                        1.0],
                                       ['V3+',
                                        'Ti4+',
                                        5.312098771970144e-05,
                                        1.0]],
                                 'sorted': [['Na1+',
                                             'Ti4+',
                                             0.00010060400812977031,
                                             1.0],
                                            ['Zn2+',
                                             'Ti4+',
                                             8.56373996146833e-05,
                                             1.0]]}}

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

        CM = self.cation_mutator

        poss_n_type_cat, poss_p_type_cat = self._get_dopants(cations, "cation")
        poss_n_type_an, poss_p_type_an = self._get_dopants(anions, "anion")

        n_type_cat, p_type_cat, n_type_an, p_type_an = [], [], [], []
        for cation in cations:
            cation_charge = utilities.parse_spec(cation)[1]

            for n_specie in poss_n_type_cat:
                n_specie_charge = utilities.parse_spec(n_specie)[1]
                if cation_charge >= n_specie_charge:
                    continue
                if CM.sub_prob(cation, n_specie) > self.threshold:
                    n_type_cat.append(
                        [
                            n_specie,
                            cation,
                            self._calculate_species_sim_prob(cation, n_specie),
                        ]
                    )

            for p_specie in poss_p_type_cat:
                p_specie_charge = utilities.parse_spec(p_specie)[1]
                if cation_charge <= p_specie_charge:
                    continue
                if CM.sub_prob(cation, p_specie) > self.threshold:
                    p_type_cat.append(
                        [
                            p_specie,
                            cation,
                            self._calculate_species_sim_prob(cation, p_specie),
                        ]
                    )
        for anion in anions:
            anion_charge = utilities.parse_spec(anion)[1]

            for n_specie in poss_n_type_an:
                n_specie_charge = utilities.parse_spec(n_specie)[1]
                if anion_charge >= n_specie_charge:
                    continue
                if CM.sub_prob(anion, n_specie) > self.threshold:
                    n_type_an.append(
                        [
                            n_specie,
                            anion,
                            self._calculate_species_sim_prob(anion, n_specie),
                        ]
                    )
            for p_specie in poss_p_type_an:
                p_specie_charge = utilities.parse_spec(p_specie)[1]
                if anion_charge <= p_specie_charge:
                    continue
                if CM.sub_prob(anion, p_specie) > self.threshold:
                    p_type_an.append(
                        [
                            p_specie,
                            anion,
                            self._calculate_species_sim_prob(anion, p_specie),
                        ]
                    )
        dopants_lists = [n_type_cat, p_type_cat, n_type_an, p_type_an]

        # sort by probability
        for dopants_list in dopants_lists:
            dopants_list.sort(key=lambda x: x[2], reverse=True)

        self.len_list = 3
        if get_selectivity:
            self.len_list = 4
            for i in range(len(dopants_lists)):
                sub = "cation"
                if i > 1:
                    sub = "anion"
                dopants_lists[i] = self._get_selectivity(
                    dopants_lists[i], cations, sub
                )

        # if groupby
        groupby_lists = [
            dict()
        ] * 4  # create list of empty dict length of 4 (n-cat, p-cat, n-an, p-an)
        # in case group_by_charge = False
        if group_by_charge:
            for i, dl in enumerate(dopants_lists):
                # groupby first element charge
                dl = sorted(dl, key=lambda x: utilities.parse_spec(x[0])[1])
                grouped_data = groupby(
                    dl, key=lambda x: utilities.parse_spec(x[0])[1]
                )
                grouped_top_data = {
                    str(k): list(g)[:num_dopants] for k, g in grouped_data
                }
                groupby_lists[i] = grouped_top_data
                del grouped_data

        # select top n elements
        dopants_lists = [
            dopants_list[:num_dopants] for dopants_list in dopants_lists
        ]

        keys = [
            "n-type cation substitutions",
            "p-type cation substitutions",
            "n-type anion substitutions",
            "p-type anion substitutions",
        ]

        self.results = self._merge_dicts(keys, dopants_lists, groupby_lists)

        # return the top (num_dopants) results for each case
        return self.results

    def plot_dopants(self, cmap: str = "YlOrRd") -> None:
        """
        Plot the dopant suggestions using the periodic table heatmap.
        Args:
            cmap (str): The colormap to use for the heatmap.
        Returns:
            None
        """
        assert (
            self.results
        ), "Dopants are not calculated. Run get_dopants first."

        for dopant_type, dopants in self.results.items():
            # due to selectivity option
            if self.len_list == 3:
                dict_results = {
                    utilities.parse_spec(x)[0]: y
                    for x, _, y in dopants.get("sorted")
                }
            else:
                dict_results = {
                    utilities.parse_spec(x)[0]: y
                    for x, _, y, _ in dopants.get("sorted")
                }
            plotting.periodic_table_heatmap(
                elemental_data=dict_results,
                cmap=cmap,
                blank_color="gainsboro",
                edge_color="white",
            )

    def _format_number(self, num_str):
        num = int(num_str)
        sign = "+" if num >= 0 else "-"
        return f"{abs(num)}{sign}"

    @property
    def to_table(self):
        if not self.results:
            print("No data available")
            return
        if self.use_probability:
            headers = ["Rank", "Dopant", "Host", "Probability", "Selectivity"]
        else:
            headers = ["Rank", "Dopant", "Host", "Similarity", "Selectivity"]
        for dopant_type, dopants in self.results.items():
            print("\033[91m" + str(dopant_type) + "\033[0m")
            for k, v in dopants.items():
                kind = k if k == "sorted" else self._format_number(k)
                print("\033[96m" + str(kind) + "\033[0m")
                enumerated_data = [
                    [i + 1] + sublist for i, sublist in enumerate(v)
                ]
                print(
                    tabulate(
                        enumerated_data,
                        headers=headers[: self.len_list + 1],
                        tablefmt="grid",
                    ),
                    end="\n\n",
                )
