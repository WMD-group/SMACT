"""Utility functions for handling oxidation states."""

from __future__ import annotations

from os import path

import pandas as pd

from smact import Element, data_directory, ordered_elements
from smact.utils.species import unparse_spec


class ICSD24OxStatesFilter:
    """Class to handle filtering the ICSD 24 oxidation states list.

    The ICSD 24 oxidation states list is a list of oxidation states for each element in the ICSD 24 database.

    Attributes:
        ox_states_df (pd.DataFrame): The ICSD 24 oxidation states list as a DataFrame.
    """

    def __init__(self):
        """Initialise the ICSD 24 oxidation states list."""
        self.ox_states_df = pd.read_json(path.join(data_directory, "oxidation_states_icsd24_counts.json"))

    def filter(
        self,
        consensus: int = 3,
        include_zero: bool = False,
        commonality: str | float = "low",
    ):
        """Filter the ICSD 24 oxidation states list by a threshold.

        Args:
            consensus (int): Minimum number of occurrences in literature for an ion to be considered valid. Default is 3.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.
            commonality (str): Excludes species below a certain proportion of appearances in literature with respect to the total number of reports of a given element (after the consensus threshold has been applied). "low" includes all species, "medium" excludes rare species below 10% occurrence, and "high" excludes non-majority species below 50% occurrence. "main" selects the species with the highest occurrence for a given element. Users may also specify their own threshold (float or int). Default is "low".

        Returns:
            pd.DataFrame: The filtered oxidation states list as a DataFrame.
        """
        commonality_map = {"low": 0, "medium": 10, "high": 50}
        commonality_threshold = 0

        if isinstance(commonality, str):
            if commonality == "main":
                pass
            else:
                commonality_threshold = commonality_map.get(commonality)
        elif isinstance(commonality, int | float):
            commonality_threshold = commonality
        else:
            raise TypeError("commonality must be a string ('low', 'medium', 'high', 'main'), a float or an integer")

        if not include_zero:
            filtered_df = self.ox_states_df[self.ox_states_df["oxidation_state"] != 0].reset_index(drop=True)
        else:
            filtered_df = self.ox_states_df

        filtered_df = filtered_df[
            (filtered_df["results_count"] >= consensus) & (filtered_df["results_count"] != 0)
        ].reset_index(drop=True)

        element_totals = filtered_df.groupby("element")["results_count"].transform("sum")
        filtered_df["species_proportion (%)"] = filtered_df["results_count"] / element_totals * 100

        if commonality == "main":
            filtered_df = (
                filtered_df.loc[filtered_df.groupby("element")["species_proportion (%)"].idxmax()]
            ).reset_index(drop=True)
        else:
            filtered_df = filtered_df[filtered_df["species_proportion (%)"] >= commonality_threshold].reset_index(
                drop=True
            )

        summary_df = (
            filtered_df.groupby("element").apply(self._filter_oxidation_states, 0, include_groups=False).reset_index()
        )
        summary_df.columns = ["element", "oxidation_state"]
        summary_df["atomic_number"] = summary_df["element"].apply(lambda x: Element(x).number)
        return summary_df.sort_values("atomic_number").drop(columns="atomic_number").reset_index(drop=True)

    def get_species_list(
        self,
        consensus: int = 3,
        include_zero: bool = False,
        include_one_oxidation_state: bool = False,
        commonality: str = "low",
    ):
        """Get the filtered ICSD 24 oxidation states list as a list of species.

        Args:
            consensus (int): Minimum number of occurrences in literature for an ion to be considered valid. Default is 3.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.
            include_one_oxidation_state (bool): Include oxidation states +1 and -1 in the filtered list or include as + and - signs. Default is False.
            commonality (str): Excludes species below a certain proportion of appearances in literature with respect to the total number of reports of a given element (after the consensus threshold has been applied). "low" includes all species, "medium" excludes rare species below 10% occurrence, and "high" excludes non-majority species below 50% occurrence. "main" selects the species with the highest occurrence for a given element. Users may also specify their own threshold (float or int). Default is "low".

        Returns:
            list: The filtered oxidation states list as a list of species.
        """
        filtered_df = self.filter(consensus, include_zero, commonality)
        species_list = []
        for _, row in filtered_df.iterrows():
            ox_states = row["oxidation_state"].split(" ")
            for ox_state in ox_states:
                try:
                    species_list.append(
                        unparse_spec(
                            (row["element"], int(ox_state)),
                            include_one=include_one_oxidation_state,
                        )
                    )
                except ValueError:
                    continue
        return species_list

    def get_species_occurrences_df(
        self,
        consensus: int = 3,
        include_one_oxidation_state: bool = False,
        sort_by_occurrences: bool = True,
        include_zero: bool = False,
    ):
        """Get the ICSD 24 oxidation states list as a dataframe of species with their occurrences.

        Args:
            consensus (int): Minimum number of occurrences in literature for an ion to be considered valid. Default is 3.
            include_one_oxidation_state (bool): Include oxidation states +1 and -1 in the species or include as + and - signs. Default is False.
            sort_by_occurrences (bool): Sort the species list by occurrences. Default is True.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.

        Returns:
            dataframe: The species list as a dataframe of species with their occurrences.
        """
        if not include_zero:
            species_occurrences_df = self.ox_states_df[self.ox_states_df["oxidation_state"] != 0].reset_index(drop=True)
        else:
            species_occurrences_df = self.ox_states_df

        species_occurrences_df = species_occurrences_df[
            (species_occurrences_df.results_count >= consensus)
        ].reset_index(drop=True)
        species_occurrences_df["species"] = species_occurrences_df.apply(
            lambda x: unparse_spec(
                (x["element"], x["oxidation_state"]),
                include_one=include_one_oxidation_state,
            ),
            axis=1,
        )
        species_occurrences_df = species_occurrences_df[["element", "species", "results_count"]]
        element_totals = species_occurrences_df.groupby("element")["results_count"].transform("sum")
        species_occurrences_df["species_proportion (%)"] = (
            species_occurrences_df["results_count"] / element_totals * 100
        )
        if sort_by_occurrences:
            return species_occurrences_df.sort_values("results_count", ascending=False).reset_index(drop=True)
        return species_occurrences_df

    def write(
        self,
        filename: str,
        comment: str | None = None,
        consensus: int = 3,
        include_zero: bool = False,
        commonality: str = "low",
    ):
        """Write the filtered ICSD 24 oxidation states list to a SMACT-compatible oxidation states txt file.

        Args:
            filename (str): The filename to write the filtered oxidation states list to.
            comment (str): A comment to include in the txt file. Default is None.
            consensus (int): Minimum number of occurrences in literature for an ion to be considered valid. Default is 3.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.
            commonality (str): Excludes species below a certain proportion of appearances in literature with respect to the total number of reports of a given element (after the consensus threshold has been applied). "low" includes all species, "medium" excludes rare species below 10% occurrence, and "high" excludes non-majority species below 50% occurrence. "main" selects the species with the highest occurrence for a given element. Users may also specify their own threshold (float or int). Default is "low".
        """
        filtered_df = self.filter(consensus, include_zero, commonality)
        # Convert the DataFrame to the require format
        summary_dict = filtered_df.set_index("element")["oxidation_state"].to_dict()
        all_elements = ordered_elements(1, 103)
        final_summary = []

        for element in all_elements:
            oxidation_states = summary_dict.get(element, "")

            if oxidation_states:
                final_summary.append(f"{element} {oxidation_states}".strip())
            else:
                final_summary.append(element)
        # Write the filtered oxidation states list to a txt file
        if not filename.endswith(".txt"):
            filename += ".txt"
        with open(filename, "w") as f:
            f.write(
                f"#\n# Oxidation state set\n# Source: ICSD (2024), filtered for {commonality} commonality of reports\n#\n"
            )
            if comment:
                f.write(f"# {comment}\n#\n")
            if include_zero:
                f.write("# Includes oxidation state 0\n#\n")
            for line in final_summary:
                f.write(line + "\n")

    def _filter_oxidation_states(self, group: pd.DataFrame.groupby, threshold: int):
        """Filter the oxidation states list by a threshold."""
        filtered_states = group[group["results_count"] >= threshold]

        return " ".join(map(str, sorted(filtered_states["oxidation_state"])))
