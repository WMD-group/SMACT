"""The dopant prediction module facilitates high-throughput prediction of p-type and n-type dopants
in multi-component solids. The search and ranking process is based on electronic filters
(e.g. accessible oxidation states) and chemical filters (e.g. difference in ionic radius).
"""

from __future__ import annotations

import logging
import os
import warnings
from itertools import groupby

import numpy as np
from pymatgen.util import plotting
from tabulate import tabulate

from smact import data_directory
from smact.structure_prediction import mutation, utilities

logger = logging.getLogger(__name__)

__all__ = [
    "SKIPSPECIES_COSINE_SIM_PATH",
    "SPECIES_M3GNET_MP2023_EFORM_COSINE_PATH",
    "SPECIES_M3GNET_MP2023_GAP_COSINE_PATH",
    "Doper",
]

_DEFAULT_LAMBDA_THRESHOLD = -5.0
_SELECTIVITY_WEIGHT = 0.25
_SIMILARITY_WEIGHT = 1 - _SELECTIVITY_WEIGHT
_PROBABILITY_INDEX = 2
_SIMILARITY_INDEX = 3
_SELECTIVITY_INDEX = 4
_COMBINED_SCORE_INDEX = 5
_NUM_DOPANT_TYPES = 4

SKIPSPECIES_COSINE_SIM_PATH = os.path.join(
    data_directory,
    "species_rep/skipspecies_20221028_319ion_dim200_cosine_similarity.json",
)
SPECIES_M3GNET_MP2023_EFORM_COSINE_PATH = os.path.join(
    data_directory, "species_rep/ion_embedding_M3GNet-MP-2023.11.1-oxi-Eform_cosine_similarity.json"
)

SPECIES_M3GNET_MP2023_GAP_COSINE_PATH = os.path.join(
    data_directory, "species_rep/ion_embedding_M3GNet-MP-2023.11.1-oxi-band_gap_cosine_similarity.json"
)

# Backward-compatible alias — the original name had a double-S typo.
SKIPSSPECIES_COSINE_SIM_PATH = SKIPSPECIES_COSINE_SIM_PATH


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
            self.cation_mutator = mutation.CationMutator.from_json(SKIPSPECIES_COSINE_SIM_PATH)
        elif embedding == "M3GNet-MP-2023.11.1-oxi-Eform":
            self.cation_mutator = mutation.CationMutator.from_json(SPECIES_M3GNET_MP2023_EFORM_COSINE_PATH)
        elif embedding == "M3GNet-MP-2023.11.1-oxi-band_gap":
            self.cation_mutator = mutation.CationMutator.from_json(SPECIES_M3GNET_MP2023_GAP_COSINE_PATH)
        elif filepath:
            self.cation_mutator = mutation.CationMutator.from_json(filepath)
        else:
            # Default to Hautier data-mined lambda values from pymatgen
            self.cation_mutator = mutation.CationMutator.from_json()
        self.possible_species = list(self.cation_mutator.specs)
        if self.cation_mutator.alpha is not None:
            self.lambda_threshold = self.cation_mutator.alpha("X", "Y")
            self.threshold = 1 / self.cation_mutator.Z * np.exp(self.cation_mutator.alpha("X", "Y"))
        else:
            self.lambda_threshold = _DEFAULT_LAMBDA_THRESHOLD
            self.threshold = 1 / self.cation_mutator.Z * np.exp(_DEFAULT_LAMBDA_THRESHOLD)
        self.use_probability = use_probability
        self.results = None

    def _get_selectivity(
        self,
        data_list: list[list],
        cations: list[str],
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
            if len(dopants) != _SELECTIVITY_INDEX + 1:  # pragma: no cover
                raise RuntimeError(
                    f"Dopant list has unexpected length {len(dopants)} (expected {_SELECTIVITY_INDEX + 1}). "
                    "This is an internal error; please report it."
                )
        return data

    def _merge_dicts(self, keys, dopants_list, groupby_list, sort_idx: int = 2):
        merged_dict = dict()
        for k, dopants, group in zip(keys, dopants_list, groupby_list, strict=False):
            merged_values = dict()
            merged_values["sorted"] = dopants
            for key, value in group.items():
                merged_values[key] = sorted(value, key=lambda x: x[sort_idx], reverse=True)
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

    def _collect_dopants(
        self,
        host_ions: list[str],
        candidates: list[str],
        cation_mutator: mutation.CationMutator,
        charge_comparison: str,
    ) -> list:
        """Collect dopants that pass the threshold filter.

        Args:
            host_ions: Host ions to substitute.
            candidates: Candidate dopant species.
            cation_mutator: The CationMutator instance for scoring.
            charge_comparison: "n_type" (candidate charge > host) or "p_type" (candidate charge < host).

        Returns:
            List of [dopant_species, host_ion, probability, lambda] entries.
        """
        results = []
        for host in host_ions:
            host_charge = utilities.parse_spec(host)[1]
            for candidate in candidates:
                candidate_charge = utilities.parse_spec(candidate)[1]
                if charge_comparison == "n_type" and host_charge >= candidate_charge:
                    continue
                if charge_comparison == "p_type" and host_charge <= candidate_charge:
                    continue
                prob = cation_mutator.sub_prob(host, candidate)
                lam = cation_mutator.get_lambda(host, candidate)
                if self.use_probability:
                    if prob <= self.threshold:
                        continue
                elif lam <= self.lambda_threshold:
                    continue
                results.append([candidate, host, prob, lam])
        return results

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
            except (AttributeError, ValueError):
                warnings.warn(f"Could not parse charge for ion '{ion}'; skipping.", stacklevel=2)

        poss_n_type_cat, poss_p_type_cat = self._get_dopants(cations, "cation")
        poss_n_type_an, poss_p_type_an = self._get_dopants(anions, "anion")

        n_type_cat = self._collect_dopants(cations, poss_n_type_cat, self.cation_mutator, "n_type")
        p_type_cat = self._collect_dopants(cations, poss_p_type_cat, self.cation_mutator, "p_type")
        n_type_an = self._collect_dopants(anions, poss_n_type_an, self.cation_mutator, "n_type")
        p_type_an = self._collect_dopants(anions, poss_p_type_an, self.cation_mutator, "p_type")
        dopants_lists = [n_type_cat, p_type_cat, n_type_an, p_type_an]

        # sort by probability or lambda depending on use_probability flag
        sort_idx = _PROBABILITY_INDEX if self.use_probability else _SIMILARITY_INDEX
        for dopants_list in dopants_lists:
            dopants_list.sort(key=lambda x: x[sort_idx], reverse=True)

        self.len_list = _NUM_DOPANT_TYPES
        if get_selectivity:
            self.len_list = _NUM_DOPANT_TYPES + 2
            for i in range(len(dopants_lists)):
                sub = "cation"
                if i > 1:
                    sub = "anion"
                dopants_lists[i] = self._get_selectivity(dopants_lists[i], cations, sub)

            for dopants_list in dopants_lists:
                for dopant in dopants_list:
                    similarity = dopant[_SIMILARITY_INDEX]
                    selectivity = dopant[_SELECTIVITY_INDEX]
                    combined_score = self._calculate_combined_score(similarity, selectivity)
                    dopant.append(combined_score)

            # sort by combined score
            for dopants_list in dopants_lists:
                dopants_list.sort(key=lambda x: x[_COMBINED_SCORE_INDEX], reverse=True)

        # if groupby
        groupby_lists = [dict() for _ in range(_NUM_DOPANT_TYPES)]  # create list of empty dict (n-cat, p-cat, n-an, p-an)
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

        # When selectivity is computed the final sort is by combined score (idx 5);
        # otherwise use the same probability-vs-lambda index as the outer sort.
        effective_sort_idx = _COMBINED_SCORE_INDEX if get_selectivity else sort_idx
        self.results = self._merge_dicts(keys, dopants_lists, groupby_lists, effective_sort_idx)

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
        if not self.results:
            raise RuntimeError("Dopants are not calculated. Run get_dopants first.")

        _plot_value_index = {
            "probability": _PROBABILITY_INDEX,
            "similarity": _SIMILARITY_INDEX,
            "selectivity": _SELECTIVITY_INDEX,
        }
        for dopants in self.results.values():
            idx = _plot_value_index.get(plot_value, _COMBINED_SCORE_INDEX)
            dict_results = {utilities.parse_spec(row[0])[0]: row[idx] for row in dopants.get("sorted")}
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
        return _SIMILARITY_WEIGHT * similarity + _SELECTIVITY_WEIGHT * selectivity

    @property
    def to_table(self) -> str:
        """
        Format the dopant suggestions as a table string.

        Returns:
        -------
            str: Formatted table of dopant suggestions.

        """
        if not self.results:
            logger.warning("No data available")
            return ""
        headers = ["Rank", "Dopant", "Host", "Probability", "Similarity", "Selectivity", "Combined"]

        parts: list[str] = []
        for dopant_type, dopants in self.results.items():
            parts.append("\033[91m" + str(dopant_type) + "\033[0m")
            for k, v in dopants.items():
                kind = k if k == "sorted" else self._format_number(k)
                parts.append("\033[96m" + str(kind) + "\033[0m")
                enumerated_data = [[i + 1, *sublist] for i, sublist in enumerate(v)]
                parts.append(
                    tabulate(
                        enumerated_data,
                        headers=headers[: self.len_list + 1],
                        tablefmt="grid",
                    )
                )
                parts.append("")
        return "\n".join(parts)
