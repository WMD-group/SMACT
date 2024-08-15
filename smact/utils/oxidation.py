"""Utility functions for handling oxidation states."""
import json
import os
from os import path

import pandas as pd

# get correct path for datafiles when called from another directory
from smact import data_directory, data_loader, ordered_elements


class ICSD24OxStatesFilter:
    """Class to handle filtering the ICSD 24 oxidation states list."""

    def __init__(self):
        self.ox_states_df = pd.read_json(
            path.join(data_directory, "oxidation_states_icsd24_counts.json")
        )

    def filter(self, threshold: int, include_zero: bool = False):
        """Filter the ICSD 24 oxidation states list by a threshold.

        Args:
            threshold (int): The threshold for filtering the oxidation states list.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.

        Returns:
            pd.DataFrame: The filtered oxidation states list as a DataFrame.
        """
        if not include_zero:
            filtered_df = self.ox_states_df[
                self.ox_states_df["oxidation_state"] != 0
            ].reset_index(drop=True)
        else:
            filtered_df = self.ox_states_df
        summary_df = (
            filtered_df[(filtered_df["results_count"] > 0)]
            .groupby("element")
            .apply(self._filter_oxidation_states, threshold)
            .reset_index()
        )
        summary_df.columns = ["element", "oxidation_state"]
        return summary_df

    def write(
        self,
        filename: str | os.PathLike,
        threshold: int,
        include_zero: bool = False,
        comment: str = None,
    ):
        """Write the filtered ICSD 24 oxidation states list to a SMACT-compatible oxidation states txt file.

        Args:
            filename (str | os.PathLike): The filename to write the filtered oxidation states list to.
            threshold (int): The threshold for filtering the oxidation states list.
            include_zero (bool): Include oxidation state of zero in the filtered list. Default is False.
            comment (str): A comment to include in the txt file. Default is None.
        """
        filtered_df = self.filter(threshold, include_zero)
        # Convert the DataFrame to the require format
        summary_dict = filtered_df.set_index("element")[
            "oxidation_state"
        ].to_dict()
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
                f"#\n# Oxidation state set\n# Source: ICSD (2024), filtered for > {threshold} reports\n#\n"
            )
            if comment:
                f.write(f"# {comment}\n")
            for line in final_summary:
                f.write(line + "\n")

    def _filter_oxidation_states(
        self, group: pd.DataFrame.groupby, threshold: int
    ):
        """Filter the oxidation states list by a threshold."""
        filtered_states = group[group["results_count"] > threshold]

        return " ".join(map(str, sorted(filtered_states["oxidation_state"])))
