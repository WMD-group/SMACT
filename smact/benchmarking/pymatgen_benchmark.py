"""Benchmarking functions for pymatgen."""

from itertools import combinations_with_replacement as cwr

from pymatgen.analysis.structure_prediction.substitution_probability import (
    SubstitutionProbability,
)

from .utilities import timeit


class ProbabilityBenchmarker:
    """Benchmarking tests for pymatgen SubstitutionProbability."""

    @timeit
    def run_tests(self):
        """Run all tests."""
        self.__sp_setup()
        self.__pair_corr()

    @timeit
    def __sp_setup(self):
        """Set up SubstitutionProbability."""
        self.sp = SubstitutionProbability()

    @timeit
    def __pair_corr(self):
        """Get pair correlation."""
        pairs = cwr(self.sp.species, 2)

        for s1, s2 in pairs:
            self.sp.pair_corr(s1, s2)


@timeit(delim=True, n=100)
def probability_test_run():
    ProbabilityBenchmarker().run_tests()
