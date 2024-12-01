"""Utility functions for handling oxidation states."""

from __future__ import annotations

from os import path

import pandas as pd

from smact import Element, data_directory, ordered_elements
from smact.structure_prediction.utilities import unparse_spec


class ICSD24OxStatesFilter:
    """Class to handle filtering the ICSD 24 oxidation states list.

    The ICSD 24 oxidation states list is a list of oxidation states for each element in the ICSD 24 database.

    Attributes:
        ox_states_df (pd.DataFrame): The ICSD 24 oxidation states list as a DataFrame.
    """

    def __init__(self):
        """Initialise the ICSD 24 oxidation states list."""
        self.ox_states_df = pd.read_json(path.join(data_directory, "oxidation_states_icsd24_counts.json"))

    def filter(self, threshold: int, include_zero: bool = False):
        """Filter the ICSD 24 oxidation states list by a threshold.

        Args:
            threshold (int): The threshold for filtering the oxidation states list.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.

        Returns:
            pd.DataFrame: The filtered oxidation states list as a DataFrame.
        """
        if not include_zero:
            filtered_df = self.ox_states_df[self.ox_states_df["oxidation_state"] != 0].reset_index(drop=True)
        else:
            filtered_df = self.ox_states_df
        summary_df = (
            filtered_df[(filtered_df["results_count"] > 0)]
            .groupby("element")
            .apply(self._filter_oxidation_states, threshold, include_groups=False)
            .reset_index()
        )
        summary_df.columns = ["element", "oxidation_state"]
        # Sort the elements by atomic number
        summary_df["atomic_number"] = summary_df["element"].apply(lambda x: Element(x).number)
        return summary_df.sort_values("atomic_number").drop(columns="atomic_number").reset_index(drop=True)

    def get_species_list(
        self,
        threshold: int,
        include_zero: bool = False,
        include_one_oxidation_state: bool = False,
    ):
        """Get the filtered ICSD 24 oxidation states list as a list of species.

        Args:
            threshold (int): The threshold for filtering the oxidation states list.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.
            include_one_oxidation_state (bool): Include oxidation states +1 and -1 in the filtered list or include as + and - signs. Default is False.

        Returns:
            list: The filtered oxidation states list as a list of species.
        """
        filtered_df = self.filter(threshold, include_zero)
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
        include_one_oxidation_state: bool = False,
        sort_by_occurrences: bool = True,
    ):
        """Get the ICSD 24 oxidation states list as a dataframe of species with their occurrences.

        Args:
            include_one_oxidation_state (bool): Include oxidation states +1 and -1 in the species or include as + and - signs. Default is False.
            sort_by_occurrences (bool): Sort the species list by occurrences. Default is True.

        Returns:
            dataframe: The species list as a dataframe of species with their occurrences.
        """
        species_occurrences_df = self.ox_states_df[
            (self.ox_states_df.results_count > 0) & (self.ox_states_df.oxidation_state != 0)
        ].reset_index(drop=True)
        species_occurrences_df["species"] = species_occurrences_df.apply(
            lambda x: unparse_spec(
                (x["element"], x["oxidation_state"]),
                include_one=include_one_oxidation_state,
            ),
            axis=1,
        )
        species_occurrences_df = species_occurrences_df[["species", "results_count"]]
        if sort_by_occurrences:
            return species_occurrences_df.sort_values("results_count", ascending=False).reset_index(drop=True)
        return species_occurrences_df

    def write(
        self,
        filename: str,
        threshold: int,
        include_zero: bool = False,
        comment: str | None = None,
    ):
        """Write the filtered ICSD 24 oxidation states list to a SMACT-compatible oxidation states txt file.

        Args:
            filename (str): The filename to write the filtered oxidation states list to.
            threshold (int): The threshold for filtering the oxidation states list.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.
            comment (str): A comment to include in the txt file. Default is None.
        """
        filtered_df = self.filter(threshold, include_zero)
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
            f.write(f"#\n# Oxidation state set\n# Source: ICSD (2024), filtered for > {threshold} reports\n#\n")
            if comment:
                f.write(f"# {comment}\n#\n")
            if include_zero:
                f.write("# Includes oxidation state 0\n#\n")
            for line in final_summary:
                f.write(line + "\n")

    def _filter_oxidation_states(self, group: pd.DataFrame.groupby, threshold: int):
        """Filter the oxidation states list by a threshold."""
        filtered_states = group[group["results_count"] > threshold]

        return " ".join(map(str, sorted(filtered_states["oxidation_state"])))
