"""Probability models for species substitution.

Implements base class :class:`SubstitutionModel`,
which can be extended to allow for development of new
lambda tables. An example of such an extension,
:class:`RadiusModel`, is also implemented.

Todo:
    * Allow for parallelism in lambda table calculations
      by implementing a `sub_probs` abstractmethod
      that :meth:`SubstitutionModel.gen_lambda` uses,
      if available.

"""

import abc
import os
from itertools import combinations_with_replacement
from typing import Dict, List, Optional

import pandas as pd

from smact import data_directory

from .utilities import parse_spec


class SubstitutionModel(abc.ABC):
    """Abstract base class for substitution models."""

    @abc.abstractmethod
    def sub_prob(self, s1: str, s2: str) -> float:
        """Calculate the probability of substituting species s1 for s2.

        Args:
            s1: The species being substituted.
            s2: The species substituting.

        Returns:
            The probability of substitution.

        """

    def gen_lambda(self, species: List[str]) -> pd.DataFrame:
        """Generate a lambda table for a list of species.

        Args:
            species: A list of species strings.

        Returns:
            A pivot table-style DataFrame containing lambda values
            for every possible species pair.

        """
        pairs = combinations_with_replacement(species, 2)

        lambda_tab = []
        for s1, s2 in pairs:
            prob = self.sub_prob(s1, s2)
            lambda_tab.append((s1, s2, prob))
            if s1 != s2:
                lambda_tab.append((s2, s1, prob))

        df = pd.DataFrame(lambda_tab)
        return df.pivot(index=0, columns=1, values=2)


class RadiusModel(SubstitutionModel):
    """Substitution probability model based on Shannon radii."""

    def __init__(self):
        r"""Parse Shannon radii data file.

        Also calculates "spring constant", _k_, based on maximum
        difference in Shannon radii:

        .. math::
            k = \Delta r_\mathrm{max}^{-2}.

        """
        shannon_file = os.path.join(data_directory, "shannon_radii.csv")

        self.shannon_data = pd.read_csv(shannon_file, index_col=0)

        self.k = (
            self.shannon_data["ionic_radius"].max()
            - self.shannon_data["ionic_radius"].min()
        ) ** -2

    def sub_prob(self, s1, s2):
        r"""Calculate the probability of substituting species s1 for s2.

        Based on the difference in Shannon radii, the probability is
        assumed to be:

        .. math::
            p = 1 - k \Delta r^2.

        Args:
            s1: The species being substituted.
            s2: The species substituting.

        Returns:
            The probability of substitution.

        """
        spec1 = parse_spec(s1)
        spec2 = parse_spec(s2)

        try:
            ele1_rows = self.shannon_data.loc[spec1[0]]
            ele2_rows = self.shannon_data.loc[spec2[0]]
        except KeyError as e:
            raise KeyError(f"Element not in Shannon radius data file: {e}")

        spec1_rows = ele1_rows.loc[ele1_rows["charge"] == spec1[1]]
        spec2_rows = ele2_rows.loc[ele2_rows["charge"] == spec2[1]]

        # Get mean so we don't need coordination information
        mean_spec1_r = spec1_rows["ionic_radius"].mean()
        mean_spec2_r = spec2_rows["ionic_radius"].mean()

        # Hooke's law-style probability
        return 1 - self.k * (mean_spec1_r - mean_spec2_r) ** 2
