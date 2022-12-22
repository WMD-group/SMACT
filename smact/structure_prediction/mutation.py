"""Tools for handling ion mutation."""

import itertools
import json
import os
import re
from copy import deepcopy
from operator import itemgetter
from typing import Callable, Generator, Optional, Tuple

import numpy as np
import pandas as pd
import pymatgen.analysis.structure_prediction as pymatgen_sp

from .structure import SmactStructure
from .utilities import parse_spec


class CationMutator:
    """Handles cation mutation of SmactStructures based on substitution probability.

    Based on the algorithm presented in:
        Hautier, G., Fischer, C., Ehrlacher, V., Jain, A., and Ceder, G. (2011)
        Data Mined Ionic Substitutions for the Discovery of New Compounds.
        Inorganic Chemistry, 50(2), 656-663.
        `doi:10.1021/ic102031h <https://pubs.acs.org/doi/10.1021/ic102031h>`_

    """

    def __init__(
        self,
        lambda_df: pd.DataFrame,
        alpha: Optional[Callable[[str, str], float]] = (lambda s1, s2: -5.0),
    ):
        """Assign attributes and get lambda table.

        Args:
            lambda_df: A pandas DataFrame, with column and index labels
                as species strings and lambda values as entries.
            alpha: A function to call to fill in missing lambda values.
                The function must take the two species' strings as
                arguments, and return a floating point lambda
                value.
                Defaults to a function that unconditionally returns
                -5.0.

        """
        self.lambda_tab = lambda_df

        self.specs = set(
            itertools.chain.from_iterable(
                set(getattr(self.lambda_tab, x)) for x in ["columns", "index"]
            )
        )

        self.alpha = alpha

        # Make sure table is fully populated
        self._populate_lambda()

        self.Z = np.exp(self.lambda_tab.to_numpy()).sum()

    @staticmethod
    def from_json(
        lambda_json: Optional[str] = None,
        alpha: Optional[Callable[[str, str], float]] = (lambda s1, s2: -5.0),
    ):
        """Create a CationMutator instance from a DataFrame.

        Args:
            lambda_json (str, optional): JSON-style representation of the
                lambda table. This is a list of entries, containing pairs
                and their associated lambda values.
                Each entry is a list of [species1, species2, lambda].
                If not supplied, defaults to the lambda table
                included with pymatgen.
            alpha: See :meth:`__init__`.

        Returns:
            A :class:`CationMutator` instance.

        """
        if lambda_json is not None:
            with open(lambda_json) as f:
                lambda_dat = json.load(f)
        else:
            # Get pymatgen lambda table
            py_sp_dir = os.path.dirname(pymatgen_sp.__file__)
            pymatgen_lambda = os.path.join(py_sp_dir, "data", "lambda.json")
            with open(pymatgen_lambda) as f:
                lambda_dat = json.load(f)

            # Get rid of 'D1+' values to reflect pymatgen
            # implementation
            lambda_dat = [x for x in lambda_dat if "D1+" not in x]

        # Convert lambda table to pandas DataFrame
        lambda_dat = [tuple(x) for x in lambda_dat]
        lambda_df = pd.DataFrame(lambda_dat)

        lambda_df = lambda_df.pivot(index=0, columns=1, values=2)

        return CationMutator(lambda_df, alpha)

    def _populate_lambda(self):
        """Populate lambda table.

        Ensures no values are NaN and performs alpha calculations,
        such that an entry exists for every possible species
        combination in the lambda table.
        Also ensures lambda table symmetry.

        """
        pairs = itertools.combinations_with_replacement(self.specs, 2)

        def add_alpha(s1, s2):
            """Add an alpha value to the lambda table at both coordinates."""
            a = self.alpha(s1, s2)
            self.lambda_tab.at[s1, s2] = a
            self.lambda_tab.at[s2, s1] = a

        def mirror_lambda(s1, s2):
            """Mirror the lambda value at (s2, s1) into (s1, s2)."""
            self.lambda_tab.at[s1, s2] = self.lambda_tab.at[s2, s1]

        for s1, s2 in pairs:
            try:
                if np.isnan(self.lambda_tab.at[s1, s2]):
                    try:
                        if not np.isnan(self.lambda_tab.at[s2, s1]):
                            mirror_lambda(s1, s2)
                        else:
                            add_alpha(s1, s2)
                    except KeyError:
                        add_alpha(s1, s2)
                else:
                    mirror_lambda(s2, s1)
            except KeyError:
                try:
                    if np.isnan(self.lambda_tab.at[s2, s1]):
                        add_alpha(s1, s2)
                    else:
                        mirror_lambda(s1, s2)
                except KeyError:
                    add_alpha(s1, s2)

        # Ensure symmetry
        idx = self.lambda_tab.index
        self.lambda_tab = self.lambda_tab[idx]

    def get_lambda(self, s1: str, s2: str) -> float:
        """Get lambda values corresponding to a pair of species.

        Args:
            s1 (str): One of the species.
            s2 (str): The other species.

        Returns:
            lambda (float): The lambda value, if it exists
                in the table. Otherwise, the alpha value
                for the two species.

        """
        if {s1, s2} <= self.specs:
            return self.lambda_tab.at[s1, s2]

        return self.alpha(s1, s2)

    def get_lambdas(self, species: str) -> pd.Series:
        """Get all the lambda values associated with a species.

        Args:
            species (str): The species for which to get the
                lambda values.

        Returns:
            A pandas Series, with index-labelled lambda entries.

        """
        if not {species} <= self.specs:
            raise ValueError(f"{species} not in lambda table.")

        return self.lambda_tab.loc[species]

    @staticmethod
    def _mutate_structure(
        structure: SmactStructure,
        init_species: str,
        final_species: str,
    ) -> SmactStructure:
        """Mutate a species within a SmactStructure.

        Replaces all instances of the species within the
        structure. Every site occupied by the species
        has its label changed to the new species,
        and the list of species of the structure is
        changed to reflect the mutation.
        Stoichiometry is maintained.
        Requires the species to have the same charge.

        Note:
            Creates a deepcopy of the supplied structure,
            such that the original instance is not modified.

        Args:
            init_species (str): The species within the structure
                to mutate.
            final_species (str): The species to replace the
                initial species with.

        Returns:
            A :class:`.~SmactStructure`, with the species
                mutated.

        """
        struct_buff = deepcopy(structure)
        init_spec_tup = parse_spec(init_species)
        struct_spec_tups = list(map(itemgetter(0, 1), struct_buff.species))
        spec_loc = struct_spec_tups.index(init_spec_tup)

        final_spec_tup = parse_spec(final_species)

        # Replace species tuple
        struct_buff.species[spec_loc] = (
            *final_spec_tup,
            struct_buff.species[spec_loc][2],
        )

        # Check for charge neutrality
        if sum(x[1] * x[2] for x in struct_buff.species) != 0:
            raise ValueError("New structure is not charge neutral.")

        # Sort species again
        struct_buff.species.sort(key=itemgetter(1), reverse=True)
        struct_buff.species.sort(key=itemgetter(0))

        # Replace sites
        struct_buff.sites[final_species] = struct_buff.sites.pop(init_species)
        # And sort
        species_strs = struct_buff._format_style("{ele}{charge}{sign}").split(
            " "
        )
        struct_buff.sites = {
            spec: struct_buff.sites[spec] for spec in species_strs
        }

        return struct_buff

    @staticmethod
    def _nary_mutate_structure(
        structure: SmactStructure,
        init_species: list,
        final_species: list,
    ) -> SmactStructure:
        """Perform a n-ary mutation of a SmactStructure (n>1).
        Replaces all instances of a group of species within the structure.
        Stoichiometry is maintained.
        Charge neutrality is preserved, but the species pair do not need the same charge.

        Args:
            init_species (list): A list of species within the structure to mutate.
            final_species (list): The list of species to replace the initial species with

        """
        # Determine the number of species to mutate
        n = len(init_species)

        struct_buff = deepcopy(structure)
        init_spec_tup_list = [parse_spec(i) for i in init_species]
        struct_spec_tups = list(map(itemgetter(0, 1), struct_buff.species))
        spec_loc = [
            struct_spec_tups.index(init_spec_tup_list[i]) for i in range(n)
        ]

        final_spec_tup_list = [parse_spec(i) for i in final_species]

        # Replace species tuple
        for i in range(n):
            struct_buff.species[spec_loc[i]] = (
                *final_spec_tup_list[i],
                struct_buff.species[spec_loc[i]][2],
            )

        # Check for charge neutrality
        if sum(x[1] * x[2] for x in struct_buff.species) != 0:
            raise ValueError("New structure is not charge neutral")

        # Sort species again
        struct_buff.species.sort(key=itemgetter(1), reverse=True)
        struct_buff.species.sort(key=itemgetter(0))

        # Replace sites
        for i in range(n):
            struct_buff.sites[final_species[i]] = struct_buff.sites.pop(
                init_species[i]
            )

        # And sort
        species_strs = struct_buff._format_style("{ele}{charge}{sign}").split(
            " "
        )
        struct_buff.sites = {
            spec: struct_buff.sites[spec] for spec in species_strs
        }

        return struct_buff

    def sub_prob(self, s1: str, s2: str) -> float:
        """Calculate the probability of substitution of two species."""
        return np.exp(self.get_lambda(s1, s2)) / self.Z

    def sub_probs(self, s1: str) -> pd.Series:
        """Determine the substitution probabilities of a species with others.

        Determines the probability of substitution of the species with every
        species in the lambda table.

        """
        probs = self.get_lambdas(s1)
        probs = np.exp(probs)
        probs /= self.Z
        return probs

    def complete_sub_probs(self) -> pd.DataFrame:
        """Generate a DataFrame with all the substitution probabilities."""
        return np.exp(self.lambda_tab) / self.Z

    def complete_cond_probs(self) -> pd.DataFrame:
        """Generate a DataFrame with all the conditional substitution probabilities."""
        lambda_exp = np.exp(self.lambda_tab)
        return lambda_exp / lambda_exp.sum(axis=1)

    def complete_pair_corrs(self) -> pd.DataFrame:
        """Generate a DataFrame with all the pair correlations."""
        corr = self.complete_sub_probs()

        # Sum each row (symmetry means this is the same as column sums)
        sums = corr.sum(axis=0)
        # Stack into matrix
        mat_sums = np.vstack([sums] * len(sums))
        # Make each element the product of (row_sum * col_sum)
        mat_sums *= sums[:, None]

        corr /= mat_sums

        return corr

    def same_spec_probs(self) -> pd.Series:
        """Calculate the same species substiution probabilities."""
        return (
            np.exp(
                pd.Series(
                    np.diag(self.lambda_tab),
                    index=[self.lambda_tab.index, self.lambda_tab.columns],
                )
            )
            / self.Z
        )

    def same_spec_cond_probs(self) -> pd.Series:
        """Calculate the same species conditional substiution probabilities."""
        return np.exp(self.lambda_tab.to_numpy().diagonal()) / np.exp(
            self.lambda_tab
        ).sum(axis=0)

    def pair_corr(self, s1: str, s2: str) -> float:
        """Determine the pair correlation of two ionic species."""
        corr = self.sub_prob(s1, s2)
        corr /= self.sub_probs(s1).sum()
        corr /= self.sub_probs(s2).sum()
        return corr

    def cond_sub_prob(self, s1: str, s2: str) -> float:
        """Calculate the probability of substitution of one species with another."""
        return (
            np.exp(self.get_lambda(s1, s2))
            / np.exp(self.get_lambdas(s2)).sum()
        )

    def cond_sub_probs(self, s1: str) -> pd.Series:
        """Calculate the probabilities of substitution of a given species.

        Calculates probabilities of substitution of given species with all
        others in the lambda table.

        """
        probs = self.get_lambdas(s1)
        probs = np.exp(probs)
        probs /= np.exp(self.lambda_tab).sum()
        return probs

    def unary_substitute(
        self,
        structure: SmactStructure,
        thresh: Optional[float] = 1e-5,
    ) -> Generator[Tuple[SmactStructure, float], None, None]:
        """Find all structures with 1 substitution with probability above a threshold.

        Args:
            structure: A :class:`SmactStructure` instance from which to generate compounds.
            thresh (float): The probability threshold; discard all substitutions that have
                a probability to generate a naturally-occuring compound less than this.

        Yields:
            Tuples of (:class:`SmactStructure`, probability).

        """
        for specie in structure.get_spec_strs():
            cond_probs = self.cond_sub_probs(specie)
            likely_probs = cond_probs.loc[cond_probs > thresh]

            for new_spec, prob in likely_probs.items():
                if any(
                    [
                        new_spec == specie,
                        parse_spec(specie)[1] != parse_spec(new_spec)[1],
                    ]
                ):
                    continue
                yield (
                    self._mutate_structure(structure, specie, new_spec),
                    prob,
                )
