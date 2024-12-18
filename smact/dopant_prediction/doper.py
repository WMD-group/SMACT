"""The dopant prediction module facilitates high-throughput prediction of p-type and n-type dopants
in multi-component solids. The search and ranking process is based on electronic filters
(e.g. accessible oxidation states) and chemical filters (e.g. difference in ionic radius).
"""

from __future__ import annotations

import os
from itertools import groupby

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
SPECIES_M3GNET_MP2023_EFORM_COSINE_PATH = os.path.join(
    data_directory, "species_rep/ion_embedding_M3GNet-MP-2023.11.1-oxi-Eform_cosine_similarity.json"
)

SPECIES_M3GNET_MP2023_GAP_COSINE_PATH = os.path.join(
    data_directory, "species_rep/ion_embedding_M3GNet-MP-2023.11.1-oxi-band_gap_cosine_similarity.json"
)


class Doper:
    """
    A class to search for n & p type dopants
    Methods: get_dopants, plot_dopants.
    """

    def __init__(
        self,
        original_species: tuple[str, ...],
        filepath: str | None = None,
        embedding: str | None = None,
        use_probability: bool = True,
    ):
        """
        Initialise the `Doper` class with a tuple of species.

        Args:
        ----
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
            raise ValueError("Only one of filepath or embedding should be provided")
        if embedding and embedding not in [
            "skipspecies",
            "M3GNet-MP-2023.11.1-oxi-Eform",
            "M3GNet-MP-2023.11.1-oxi-band_gap",
        ]:
            raise ValueError(f"Embedding {embedding} is not supported")
        if embedding == "skipspecies":
            self.cation_mutator = mutation.CationMutator.from_json(SKIPSSPECIES_COSINE_SIM_PATH)
        elif embedding == "M3GNet-MP-2023.11.1-oxi-Eform":
            self.cation_mutator = mutation.CationMutator.from_json(SPECIES_M3GNET_MP2023_EFORM_COSINE_PATH)
        elif embedding == "M3GNet-MP-2023.11.1-oxi-band_gap":
            self.cation_mutator = mutation.CationMutator.from_json(SPECIES_M3GNET_MP2023_GAP_COSINE_PATH)
        elif filepath:
            self.cation_mutator = mutation.CationMutator.from_json(filepath)
        else:
            # Default to Hautier data-mined lambda values
            self.cation_mutator = mutation.CationMutator.from_json(filepath)
        self.possible_species = list(self.cation_mutator.specs)
        self.lambda_threshold = self.cation_mutator.alpha("X", "Y")
        self.threshold = 1 / self.cation_mutator.Z * np.exp(self.cation_mutator.alpha("X", "Y"))
        self.use_probability = use_probability
        self.results = None

    def _get_selectivity(
        self,
        data_list: list[smact.Element],
        cations: list[smact.Element],
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
                    sum_prob += self.cation_mutator.sub_prob(cation, selected_site)

            selectivity = sub_prob / sum_prob
            selectivity = round(selectivity, 2)
            dopants.append(selectivity)
            assert len(dopants) == 5
        return data

    def _merge_dicts(self, keys, dopants_list, groupby_list):
        merged_dict = dict()
        for k, dopants, group in zip(keys, dopants_list, groupby_list, strict=False):
            merged_values = dict()
            merged_values["sorted"] = dopants
            for key, value in group.items():
                merged_values[key] = sorted(value, key=lambda x: x[2], reverse=True)
            merged_dict[k] = merged_values
        return merged_dict

    def _get_dopants(
        self,
        specie_ions: list[str],
        ion_type: str,
    ):
        """
        Get possible dopants for a given list of elements and dopants.

        Args:
        ----
            specie_ions (List[str]): List of original species (anions or cations) as strings.
            ion_type (str): Identify which species to check.

        Returns:
        -------
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

    def get_dopants(self, num_dopants: int = 5, get_selectivity=True, group_by_charge=True) -> dict:
        """
        Get the top n dopants for each case.

        Args:
        ----
            num_dopants (int): The number of dopants to return.
            get_selectivity (bool): Whether to calculate the selectivity of the dopants.
            group_by_charge (bool): Whether to group the dopants by charge.

        Returns:
        -------
            dict: A dictionary of the top n dopants for each case.

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

            for _i, n_specie in enumerate(poss_n_type_cat):
                n_specie_charge = utilities.parse_spec(n_specie)[1]
                if cation_charge >= n_specie_charge:
                    continue
                if CM.sub_prob(cation, n_specie) > self.threshold:
                    n_type_cat.append(
                        [
                            n_specie,
                            cation,
                            self.cation_mutator.sub_prob(cation, n_specie),
                            self.cation_mutator.get_lambda(cation, n_specie),
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
                            self.cation_mutator.sub_prob(cation, p_specie),
                            self.cation_mutator.get_lambda(cation, p_specie),
                        ],
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
                            self.cation_mutator.sub_prob(anion, n_specie),
                            self.cation_mutator.get_lambda(anion, n_specie),
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
                            self.cation_mutator.sub_prob(anion, p_specie),
                            self.cation_mutator.get_lambda(anion, p_specie),
                        ]
                    )
        dopants_lists = [n_type_cat, p_type_cat, n_type_an, p_type_an]

        # sort by probability
        for dopants_list in dopants_lists:
            dopants_list.sort(key=lambda x: x[2], reverse=True)

        self.len_list = 4
        if get_selectivity:
            self.len_list = 6
            for i in range(len(dopants_lists)):
                sub = "cation"
                if i > 1:
                    sub = "anion"
                dopants_lists[i] = self._get_selectivity(dopants_lists[i], cations, sub)

            for dopants_list in dopants_lists:
                for dopant in dopants_list:
                    similarity = dopant[3]
                    selectivity = dopant[4]
                    combined_score = self._calculate_combined_score(similarity, selectivity)
                    dopant.append(combined_score)

            # sort by combined score
            for dopants_list in dopants_lists:
                dopants_list.sort(key=lambda x: x[5], reverse=True)

        # if groupby
        groupby_lists = [dict()] * 4  # create list of empty dict length of 4 (n-cat, p-cat, n-an, p-an)
        # in case group_by_charge = False
        if group_by_charge:
            for i, dl in enumerate(dopants_lists):
                # groupby first element charge
                dl = sorted(dl, key=lambda x: utilities.parse_spec(x[0])[1])
                grouped_data = groupby(dl, key=lambda x: utilities.parse_spec(x[0])[1])
                grouped_top_data = {str(k): list(g)[:num_dopants] for k, g in grouped_data}
                groupby_lists[i] = grouped_top_data
                del grouped_data

        # select top n elements
        dopants_lists = [dopants_list[:num_dopants] for dopants_list in dopants_lists]

        keys = [
            "n-type cation substitutions",
            "p-type cation substitutions",
            "n-type anion substitutions",
            "p-type anion substitutions",
        ]

        self.results = self._merge_dicts(keys, dopants_lists, groupby_lists)

        # return the top (num_dopants) results for each case
        return self.results

    def plot_dopants(self, cmap: str = "YlOrRd", plot_value: str = "probability") -> None:
        """
        Plot the dopant suggestions using the periodic table heatmap.

        Args:
        ----
            cmap (str): The colormap to use for the heatmap.
            plot_value (str): The value to plot on the heatmap. Options are "probability", "similarity" or "selectivity".

        Returns:
        -------
            None

        """
        assert self.results, "Dopants are not calculated. Run get_dopants first."

        for dopants in self.results.values():
            # due to selectivity option
            if self.len_list == 3:
                dict_results = {utilities.parse_spec(x)[0]: y for x, _, y in dopants.get("sorted")}
            elif plot_value == "probability":
                dict_results = {utilities.parse_spec(x)[0]: y for x, _, y, _, _, _ in dopants.get("sorted")}
            elif plot_value == "similarity":
                dict_results = {utilities.parse_spec(x)[0]: y for x, _, _, y, _, _ in dopants.get("sorted")}
            elif plot_value == "selectivity":
                dict_results = {utilities.parse_spec(x)[0]: y for x, _, _, _, y, _ in dopants.get("sorted")}
            else:
                dict_results = {utilities.parse_spec(x)[0]: y for x, _, _, _, _, y in dopants.get("sorted")}
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

    def _calculate_combined_score(self, similarity: float, selectivity: float) -> float:
        return (1 - 0.25) * similarity + 0.25 * selectivity

    @property
    def to_table(self):
        """
        Print the dopant suggestions in a tabular format.

        Returns:
        -------
            None

        """
        if not self.results:
            print("No data available")
            return
        # if self.use_probability:
        #     headers = ["Rank", "Dopant", "Host", "Probability", "Selectivity", "Combined"]
        # else:
        #     headers = ["Rank", "Dopant", "Host", "Similarity", "Selectivity", "Combined"]

        headers = ["Rank", "Dopant", "Host", "Probability", "Similarity", "Selectivity", "Combined"]

        for dopant_type, dopants in self.results.items():
            print("\033[91m" + str(dopant_type) + "\033[0m")
            for k, v in dopants.items():
                kind = k if k == "sorted" else self._format_number(k)
                print("\033[96m" + str(kind) + "\033[0m")
                enumerated_data = [[i + 1, *sublist] for i, sublist in enumerate(v)]
                print(
                    tabulate(
                        enumerated_data,
                        headers=headers[: self.len_list + 1],
                        tablefmt="grid",
                    ),
                    end="\n\n",
                )
